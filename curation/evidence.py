#!/usr/bin/python
"""The evidence dependency graph: what each piece of evidence is derived from and how the run used it.

A node is an evidence record (``evidence.<role>``, id ``evidence:<role>`` unless given) or an ``evidence_sources``
entry.  An edge points from a node to the sources it was derived from.  Provenance (ancestry) is kept apart from
*scope*: whether a node presently acts as a constraint or validator in the assessed interval.

A validator shares lineage with the constraining evidence when, anywhere in the two ancestries:

* the same node id occurs (``same_node``);
* two differently named nodes name the same resource, an equal SHA-256 or an equal source string, with equal
  declared selectors (``same_resource``) or with a selector missing on either side
  (``same_resource_selector_unknown``);
* two nodes name one container but declare different selectors, dataset/channel/fit/calibration identity
  (``shared_container_different_selector``): recorded as a distinction, not as dependence.

Not sharing declared lineage is *not* a proof of independence; undeclared dependencies are invisible, and a
different string or hash never establishes independence.
"""

from curation import intervals
from curation.manifest import CONSTRAINING

DEPENDENT = ('same_node', 'same_resource', 'same_resource_selector_unknown')


def build_graph(manifest, excluded=(), scopes=None):
    """``{id: node}`` from the effective manifest; ``excluded`` ids (not applicable here) are left out.

    ``scopes`` maps a node id to its effective time scope ``(start, end)`` (closed, ``None`` = unbounded side) or to
    ``False`` when the node acts nowhere in the assessed interval; an id with no entry acts everywhere.
    """
    nodes = {}
    for role, record in manifest.get('evidence', {}).items():
        identifier = record.get('id', 'evidence:' + role)
        nodes[identifier] = {'id': identifier, 'kind': 'evidence', 'role': role, 'parents': list(record.get('parents', [])),
                             'usage': list(record.get('usage', [])), 'source': record.get('source', ''),
                             'sha256': record.get('sha256'), 'selector': record.get('selector')}
    for entry in manifest.get('evidence_sources', []):
        nodes[entry['id']] = {'id': entry['id'], 'kind': 'source', 'role': None, 'parents': list(entry['parents']),
                              'usage': list(entry['usage']), 'source': entry.get('source', ''),
                              'sha256': entry.get('sha256'), 'selector': entry.get('selector'),
                              'description': entry['description']}
    for identifier in excluded:
        nodes.pop(identifier, None)
    for identifier, node in nodes.items():
        node['parents'] = [p for p in node['parents'] if p in nodes]
        scope = (scopes or {}).get(identifier)
        node['scope'] = None if scope is None else (False if scope is False else list(scope))
    return nodes


def ancestry(nodes, identifier):
    """The node and every ancestor reachable through parents."""
    seen, stack = set(), [identifier]
    while stack:
        current = stack.pop()
        if current in seen or current not in nodes:
            continue
        seen.add(current)
        stack.extend(nodes[current]['parents'])
    return seen


def acts_together(nodes, a, b):
    """Whether two nodes both act in the assessed interval and their scopes meet (unscoped acts everywhere)."""
    scope_a, scope_b = nodes[a]['scope'], nodes[b]['scope']
    if scope_a is False or scope_b is False:
        return False
    if scope_a is None or scope_b is None:
        return True
    return intervals.overlap(scope_a, scope_b) or scope_a[0] == scope_b[1] or scope_b[0] == scope_a[1]


def _identities(node):
    keys = []
    if node['sha256']:
        keys.append(('sha256', node['sha256'].lower()))
    if node['source']:
        keys.append(('source', node['source']))
    return keys


def _alias_kind(a, b):
    """How two distinct nodes naming one resource relate, from their selectors."""
    sel_a, sel_b = a['selector'], b['selector']
    if sel_a is None or sel_b is None:
        return 'same_resource_selector_unknown'
    return 'same_resource' if sel_a == sel_b else 'shared_container_different_selector'


def shared_lineage(nodes, validator):
    """Relations between ``validator`` and the constraining evidence: ``[{'with', 'kind', 'via', 'identity'}]``.

    Only constraining nodes that act in the same part of the assessed interval as the validator are compared;
    ancestry itself is never pruned (a shared parent stays visible even when its use is time-disjoint).
    """
    found = []
    mine = ancestry(nodes, validator)
    for identifier, other in sorted(nodes.items()):
        if not set(other['usage']) & set(CONSTRAINING) or not acts_together(nodes, validator, identifier):
            continue
        theirs = ancestry(nodes, identifier)
        common = sorted(mine & theirs)
        if common:
            found.append({'with': identifier, 'kind': 'same_node', 'via': common, 'identity': None})
        for node_v in sorted(mine - theirs):
            for node_c in sorted(theirs - mine):
                for key in _identities(nodes[node_v]):
                    if key in _identities(nodes[node_c]):
                        found.append({'with': identifier, 'kind': _alias_kind(nodes[node_v], nodes[node_c]),
                                      'via': [node_v, node_c], 'identity': list(key)})
    seen, unique = set(), []
    for item in found:
        marker = (item['with'], item['kind'], tuple(item['via']), tuple(item['identity'] or ()))
        if marker not in seen:
            seen.add(marker)
            unique.append(item)
    return unique


def summary(nodes):
    """A JSON-ready description of the graph, with children added."""
    children = {identifier: [] for identifier in nodes}
    for identifier, node in nodes.items():
        for parent in node['parents']:
            children[parent].append(identifier)
    return {identifier: {**node, 'children': sorted(children[identifier])} for identifier, node in sorted(nodes.items())}
