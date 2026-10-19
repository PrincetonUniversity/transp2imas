import unittest,sys,xml.etree.ElementTree as ET
sys.path.insert(0,'/mnt/data/pt_latest')
import types
sys.modules['input_provenance']=types.SimpleNamespace(assignments=lambda x: [],strip_fortran_comment=lambda x:x)
from pt_solver_regions import region_occurrence_xml,region_validity_xml,numerical_settings_xml
class TestValidityNumerics(unittest.TestCase):
 def setUp(self):
  self.plan={'boundary_states':[{'start_time':None,'values':{'XIMIN_CONF':.1,'XIMAX_CONF':.75,'XIBOUND':.8,'XNBOUND':.9,'XPHIBOUND':.95}},{'start_time':.5,'values':{'XIMIN_CONF':.2,'XIMAX_CONF':.7,'XIBOUND':.85,'XNBOUND':.92,'XPHIBOUND':.98}}], 'templates':[{'path':'pt.nml','sha256':'test','blocks':[{'start_time':-1,'options':{'newton_iterates':'100','pt_residual%res_te':'1D-3','theta_implicit':'-1D0'},'regions':{'axial':{'exb':{'version':'3'}}}},{'start_time':.5,'options':{'newton_iterates':'50'},'regions':{}}]}], 'trusted_model_assignment':True,'trdat_file':'TR.DAT','trdat_present':True,'predictive_mode':3,'model_source':'explicit_PT_TEMPLATE'}
  for s in self.plan['boundary_states']:s['origins']={k:'TR.DAT' for k in s['values']}
  for x in self.plan['templates']:x['text']='&pt_controls /'
 def test_domains_are_channel_and_time_specific(self):
  root=region_validity_xml(self.plan,'confinement'); states=root.findall('state');self.assertEqual(len(states),2)
  d={z.get('quantity'):z for z in states[0].findall('domain')}
  self.assertEqual(d['temperature'].get('xi_max'),'0.75');self.assertEqual(d['density'].get('xi_max'),'0.75');self.assertEqual(d['toroidal_momentum'].get('xi_max'),'0.75')
  self.assertEqual(states[1].get('effective_from'),'0.5')
 def test_numerics_include_versions_updates_and_residuals(self):
  root=numerical_settings_xml(self.plan)
  controls=root.findall('.//control');names={c.get('name') for c in controls}
  self.assertTrue({'newton_iterates','pt_residual%res_te','theta_implicit'}.issubset(names))
  self.assertEqual(root.find('.//regional_model_control').get('raw'),'3')
  self.assertEqual(len(root.findall('.//settings')),2)
 def test_occurrence_has_direct_validity_and_numerics(self):
  item={'region':'axial','models':{'anomalous':'glf23','neoclassical':'neoch'}}
  r=ET.fromstring(region_occurrence_xml(self.plan,item,1))
  self.assertEqual(len(r.findall('./regional_validity/state')),2)
  self.assertTrue(r.findall('.//pt_solver_numerics/input_file/settings'))
if __name__=='__main__':unittest.main()
