"""Tests for stocal.experimental.mm_kinetics module"""
import unittest
from stocal.experimental.mm_kinetics import MichaelisMenten
from stocal.tests.test_transitions import TestMassAction


class TestMMKinetics(TestMassAction):
    MichaelisMenten = MichaelisMenten

    def test_init(self):
        with self.assertRaises(ValueError):
            self.MichaelisMenten({'S': 1, 'E': 1}, {'P': 1, 'E': 1}, 'E', -1, 1)
        with self.assertRaises(ValueError):
            self.MichaelisMenten({'S': 1, 'E': 1}, {'P': 1, 'E': 1}, 'E', 1, -1)
        with self.assertRaises(ValueError):
            self.MichaelisMenten({'a': 1, 'b': 1, 'd': 1}, {'c': 1, 'd': 1}, 'd', 1, 1)
        with self.assertRaises(ValueError):
            self.MichaelisMenten({'a': 1, 'b': 1}, {'a': 1, 'd': 1, 'e': 1}, 'a', 1, 1)
        with self.assertRaises(ValueError):
            self.MichaelisMenten({'a': 1, 'b': 1}, {'c': 1, 'd': 1}, 'a', 1, 1)
        with self.assertRaises(ValueError):
            self.MichaelisMenten({'a': 1, 'b': 1}, {'c': 1, 'd': 1}, 'c', 1, 1)

    def test_equality(self):
        self.assertEqual(self.MichaelisMenten({'S': 1, 'E': 1}, {'P': 1, 'E': 1}, 'E', 1, 1),
                         self.MichaelisMenten({'S': 1, 'E': 1}, {'P': 1, 'E': 1}, 'E', 1, 1))
        self.assertEqual(self.MichaelisMenten({'E': 1, 'S': 1}, {'P': 1, 'E': 1}, 'E', 1, 1),
                         self.MichaelisMenten({'S': 1, 'E': 1}, {'P': 1, 'E': 1}, 'E', 1, 1))
        self.assertNotEqual(self.MichaelisMenten({'S': 1, 'E': 1}, {'P': 1, 'E': 1}, 'E', 1, 1),
                            self.MichaelisMenten({'S': 1, 'E': 1}, {'P': 1, 'E': 1}, 'E', 2, 1))
        self.assertNotEqual(self.MichaelisMenten({'S': 1, 'E': 1}, {'P': 1, 'E': 1}, 'E', 1, 1),
                            self.MichaelisMenten({'S': 1, 'E': 1}, {'P': 1, 'E': 1}, 'E', 1, 2))

    def test_repr(self):
        self.assertEqual(self.MichaelisMenten({'S': 1, 'E': 1}, {'P': 1, 'E': 1}, 'E', 1, 1).__repr__(),
                         '<MichaelisMenten S + E --> P + E>')

    def test_propensity_equality(self):
        k_cat = 1.5
        n_E = 10
        n_S = 10
        v_max = k_cat * n_E
        k_m = 5
        a_1 = self.MichaelisMenten({'S': 1, 'E': 1}, {'P': 1, 'E': 1}, 'E', k_cat, k_m).propensity({'S': n_S, 'E': n_E})
        a_2 = (v_max * n_S) / (k_m + n_S)
        self.assertEqual(a_1, a_2)


if __name__ == '__main__':
    unittest.main()
