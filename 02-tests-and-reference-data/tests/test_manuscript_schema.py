"""Prescribed environmental hypotheses remain optional and strictly validated."""
import copy
import unittest
from chromatin_studio.schema import default_config, preset_config, validate_config

class ManuscriptSchemaTests(unittest.TestCase):
    def test_old_config_and_new_manuscript_presets_validate(self):
        for config in [default_config(), preset_config('manuscript_locus'), preset_config('manuscript_20mb')]:
            with self.subTest(name=config['name']):
                report=validate_config(config)
                self.assertTrue(report['valid'],report['errors'])
        self.assertEqual(preset_config('manuscript_locus')['chromatin']['bp_per_bead'],200)
        large=preset_config('manuscript_20mb')
        self.assertEqual(large['chromatin']['beads']*large['chromatin']['bp_per_bead'],20_000_000)
        self.assertEqual(large['chromatin']['bp_per_bead'],10_000)

    def test_zero_stiffness_is_explicit_unconfined_control(self):
        config=default_config();config['simulation']['confinement_stiffness']=0
        report=validate_config(config)
        self.assertTrue(report['valid'])
        self.assertTrue(any('Confinement force is disabled' in w for w in report['warnings']))
        for invalid in (-1,True,float('inf'),201):
            config['simulation']['confinement_stiffness']=invalid
            self.assertFalse(validate_config(config)['valid'])

    def test_response_must_refer_to_known_body_and_positive_multiplier(self):
        config=default_config()
        for responses in ({'missing':{'activation_multiplier':1}}, {'speckle':{'activation_multiplier':0}},
                          {'speckle':{'activation_multiplier':float('nan')}}, {'speckle':{'activation_multiplier':True}},
                          {'speckle':{'off_rate':3}}, {'speckle':[]}, []):
            with self.subTest(response=responses):
                config['genes'][0]['body_responses']=responses
                self.assertFalse(validate_config(config)['valid'])

    def test_neutral_or_partial_response_is_valid_without_mutation(self):
        config=default_config();config['genes'][0]['body_responses']={'speckle':{'activation_multiplier':1}}
        before=copy.deepcopy(config)
        self.assertTrue(validate_config(config)['valid'])
        self.assertEqual(config,before)

if __name__=='__main__':unittest.main()
