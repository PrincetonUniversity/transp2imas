#!/usr/bin/python
"""The evidence dependency graph: what each piece of evidence is derived from and how the run used it.

A node is an evidence record (``evidence.<role>``, id ``evidence:<role>`` unless given) or an ``evidence_sources``
entry.  An edge points from a node to the sources it was derived from.  Independence of a validator is judged on
lineage only: it shares lineage with the constraining evidence when its ancestry (itself included) meets the
ancestry of any node used as a constraint, calibration target, initial or boundary condition, or when the two
carry the same nonempty source string or the same SHA-256.  Not sharing declared lineage is *not* a proof of
independence; undeclared dependencies cannot be seen.
"""

from curation.manifest import CONSTRAINING


def build_graph(manifest, excluded=()):
    """``{id: node}`` from the effective manifest; ``excluded`` ids (not applicable here) are left out."""
    nodes = {}
    for role, record in manifest.get('evidence', {}).items():
        identifier = record.get('id', 'evidence:' + role)
        nodes[identifier] = {'id': identifier, 'kind': 'evidence', 'role': role, 'parents': list(record.get('parents', [])),
                             'usage': list(record.get('usage', [])), 'source': record.get('source', ''),
                             'sha256': record.get('sha256')}
    for entry in manifest.get('evidence_sources', []):
        nodes[entry['id']] = {'id': entry['id'], 'kind': 'source', 'role': None, 'parents': list(entry['parents']),
                              'usage': list(entry['usage']), 'source': entry.get('source', ''),
                              'sha256': entry.get('sha256'), 'description': entry['description']}
    for identifier in excluded:
        nodes.pop(identifier, None)
    for node in nodes.values():
        node['parents'] = [p for p in node['parents'] if p in nodes]
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


def shared_lineage(nodes, validator):
    """Why ``validator`` is not independent of the constraining evidence: ``[{'with', 'via'}]`` (empty if none)."""
    mine = ancestry(nodes, validator)
    node = nodes[validator]
    found = []
    for identifier, other in nodes.items():
        if not set(other['usage']) & set(CONSTRAINING):
            continue
        common = sorted(mine & ancestry(nodes, identifier))
        if common:
            found.append({'with': identifier, 'via': common})
        elif node['source'] and node['source'] == other['source']:
            found.append({'with': identifier, 'via': ['same source: ' + node['source']]})
        elif node['sha256'] and node['sha256'].lower() == (other['sha256'] or '').lower():
            found.append({'with': identifier, 'via': ['same sha256']})
    return found


def summary(nodes):
    """A JSON-ready description of the graph, with children added."""
    children = {identifier: [] for identifier in nodes}
    for identifier, node in nodes.items():
        for parent in node['parents']:
            children[parent].append(identifier)
    return {identifier: {**node, 'children': sorted(children[identifier])} for identifier, node in sorted(nodes.items())}
