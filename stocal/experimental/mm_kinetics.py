from stocal.transitions import Reaction
from stocal.structures import multiset


class MichaelisMenten(Reaction):
    """
    Implementation of Michaelis Menten kinetics(further MM) reaction

    Since MM calculations are meant for enzyme catalyzed reactions, this method expects reactants and products that
    would follow the E + S -> E + P reaction pattern.

    The propensity is calculated based on maximum velocity of the reaction(v_max) and Michaelis constant(k_m) which
    are expected to be provided by user.
    """
    def __init__(self, reactants, products, enzyme, k_cat, k_m):
        """Data verification"""
        if k_cat < 0:
            raise ValueError("K_cat must be non negative")
        if k_m < 0:
            raise ValueError("K_m must be non negative")
        if len(reactants) != 2:
            raise ValueError("Reactants are not suitable for Michaelis Menten calculations")
        if len(products) != 2:
            raise ValueError("Products are not suitable for Michaelis Menten calculations")
        if not (list(reactants.keys())[0] == enzyme or list(reactants.keys())[1] == enzyme):
            raise ValueError("Enzyme provided is not among the reactants")
        if not (list(products.keys())[0] == enzyme or list(products.keys())[1] == enzyme):
            raise ValueError("Enzyme provided is not among the products")

        super(MichaelisMenten, self).__init__(reactants, products)
        self.enzyme = enzyme
        self.k_cat = k_cat
        if list(reactants.keys())[0] == self.enzyme:
            self.substrate = list(reactants.keys())[1]
        else:
            self.substrate = list(reactants.keys())[0]
        self.k_m = k_m

    def __repr__(self):
        try:
            return '<%s %s, %g %g>' % (
                type(self).name, self, self.enzyme, self.k_cat, self.k_m
                )
        except AttributeError:
            return super(MichaelisMenten, self).__repr__()

    def propensity(self, state):
        state = state if isinstance(state, multiset) else multiset(state)
        v_max = state[self.enzyme] * self.k_cat
        return (v_max * state[self.substrate]) / (self.k_m + state[self.substrate])

    def __eq__(self, other):
        return (
            super(MichaelisMenten, self).__eq__(other) and
            isinstance(other, MichaelisMenten) and
            self.k_m == other.k_m and
            self.k_cat == other.k_cat and
            self.enzyme == other.enzyme and
            self.substrate == other.substrate
        )

    def __hash__(self):
        if not self._hash:
            self._hash = hash((
                super(MichaelisMenten, self).__hash__(),
                self.k_cat,
                self.k_m,
                self.enzyme,
                self.substrate
            ))
        return self._hash

