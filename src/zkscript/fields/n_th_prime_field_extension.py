"""Bitcoin scripts that perform arithmetic operations in F_q^n = F_q[t] / (t^n - NON_RESIDUE)."""

from tx_engine import Script

from src.zkscript.fields.fq import Fq
from src.zkscript.fields.prime_field_extension import PrimeFieldExtension
from src.zkscript.types.stack_elements import StackFiniteFieldElement
from src.zkscript.util.utility_functions import check_order
from src.zkscript.util.utility_scripts import (
    bitmask_to_boolean_list,
    bool_to_moving_function,
    move,
    nums_to_script,
    pick,
    verify_bottom_constant,
)


class NthPrimeFieldExtension(PrimeFieldExtension):
    """Construct Bitcoin scripts that perform arithmetic operations in F_q^n = F_q[t] / (t^n - NON_RESIDUE).

    Elements in F_q^n are of the form `a_0 + a_1 * t + .. a_(n_1) t^n`, where the `a_i`s are elements of F_q.

    Attributes:
        MODULUS: The characteristic of the field F_q.
        EXTENSION_DEGREE: The extension degree over the prime field, equal to 6.
        NON_RESIDUE: The non n-th non residue used to construct the extension.
        PRIME_FIELD: The Bitcoin Script implementation of the prime field F_q.
    """

    def __init__(self, q: int, extension_degree: int, non_residue: int):
        """Initialise F_q^n.

        Args:
            q (int): The characteristic of the field F_q.
            extension_degree (int): The extension degree, i.e., n.
            non_residue (int): The n-th non residue used to construct the extension.
        """
        self.MODULUS = q
        self.EXTENSION_DEGREE = extension_degree
        self.NON_RESIDUE = non_residue
        self.PRIME_FIELD = Fq(q)

    def mul(
        self,
        take_modulo: bool,
        positive_modulo: bool = True,
        check_constant: bool | None = None,
        clean_constant: bool | None = None,
        is_constant_reused: bool | None = None,
        x: StackFiniteFieldElement | None = None,
        y: StackFiniteFieldElement | None = None,
        rolling_options: int = 3,
    ) -> Script:
        r"""Multiplication in F_q^n = F_q[t] / (t^n - NON_RESIDUE).

        Stack input:
            - stack:    [q, ..., x := (x_0, .., x_(n-1)), .., y := (y_0, .., y_(n-1)), ..]
            - altstack: []

        Stack output:
            - stack:    [q, ..., x * y := (z_0, .., z_(n-1))]
            - altstack: []

        Where:
            z_k = \sum_(i+j = k mod n, i + j < n) x_i * y_j + (\sum_(i+j = k mod n, i + j >= n) x_i * y_j) * NON_RESIDUE

        Args:
            take_modulo (bool): If `True`, the result is reduced modulo `q`.
            positive_modulo (bool): If `True` the modulo of the result is taken positive. Defaults to `True`.
            check_constant (bool | None): If `True`, check if `q` is valid before proceeding. Defaults to `None`.
            clean_constant (bool | None): If `True`, remove `q` from the bottom of the stack. Defaults to `None`.
            is_constant_reused (bool | None, optional): If `True`, `q` remains as the second-to-top element on the stack
                after execution. Defaults to `None`.
            x (StackFiniteFieldElement): The position in the stack of `x`.
            y (StackFiniteFieldElement): The position in the stack of `y`.
            rolling_options (int): Bitmask detailing which of the elements `x` and `y` should be removed from the stack
                after execution. Defaults to `3` (remove everything).

        Returns:
            Script to add two elements in F_q^n.
        """
        x = x if x is not None else StackFiniteFieldElement(2 * self.EXTENSION_DEGREE - 1, False, self.EXTENSION_DEGREE)
        y = y if y is not None else StackFiniteFieldElement(self.EXTENSION_DEGREE - 1, False, self.EXTENSION_DEGREE)
        check_order([x, y])
        is_x_rolled, is_y_rolled = bitmask_to_boolean_list(rolling_options, 2)

        out = verify_bottom_constant(self.MODULUS) if check_constant else Script()

        for k in range(self.EXTENSION_DEGREE - 1, -1, -1):
            is_last_iteration = k == 0
            x_moving_function = bool_to_moving_function(is_x_rolled) if is_last_iteration else pick
            y_moving_function = bool_to_moving_function(is_y_rolled) if is_last_iteration else pick
            # Compute \sum_(i + j = k mod n, i + j >= n) x_i * y_j
            for step in range(self.EXTENSION_DEGREE - 1, k, -1):
                out += move(
                    y.shift(1 if step < self.EXTENSION_DEGREE - 1 else 0)  # Running calculation
                    .shift(
                        -(self.EXTENSION_DEGREE - 1 - step) * is_y_rolled if is_last_iteration else 0
                    )  # Shifts from rolling
                    .extract_component(step),
                    y_moving_function,
                )
                out += move(
                    x.shift(1 if step < self.EXTENSION_DEGREE - 1 else 0)  # Running calculation
                    .shift(1)  # Moving y_step in this iteration
                    .shift(
                        -is_y_rolled * (self.EXTENSION_DEGREE - step) if is_last_iteration else 0
                    )  # Shifts from rolling/picking
                    .extract_component(k + self.EXTENSION_DEGREE - step),
                    x_moving_function,
                )
                out += Script.parse_string("OP_MUL")
                out += Script.parse_string("OP_ADD" if step < self.EXTENSION_DEGREE - 1 else "")
            if k != self.EXTENSION_DEGREE - 1:
                out += nums_to_script([self.NON_RESIDUE])
                out += Script.parse_string("OP_MUL")
            # Compute \sum_(i + j = k mod n, i + j < n) x_i * y_j
            for step in range(k, -1, -1):
                out += move(
                    y.shift(1 if k != self.EXTENSION_DEGREE - 1 and step == k else 0)  # Calculation above
                    .shift(1 if step < k else 0)  # Running calculation
                    .shift(
                        -(self.EXTENSION_DEGREE - 1 - step) * is_y_rolled if is_last_iteration else 0
                    )  # Shifts from rolling
                    .extract_component(step),
                    y_moving_function,
                )
                out += move(
                    x.shift(1 if k != self.EXTENSION_DEGREE - 1 and step == k else 0)  # Calculation above
                    .shift(1 if step < k else 0)  # Running calculation
                    .shift(1)  # Moving y_step in this iteration
                    .shift(
                        -is_y_rolled * (self.EXTENSION_DEGREE - step) - is_x_rolled * (self.EXTENSION_DEGREE - 1 - step)
                        if is_last_iteration
                        else 0
                    )  # Shifts from rolling/picking
                    .extract_component(k - step),
                    x_moving_function,
                )
                out += Script.parse_string("OP_MUL")
                out += Script.parse_string(
                    "OP_ADD" if step < k or (step == k and k != self.EXTENSION_DEGREE - 1) else ""
                )
            out += Script.parse_string("OP_TOALTSTACK" if k > 0 else "")

        out += (
            self.take_modulo(
                positive_modulo=positive_modulo, clean_constant=clean_constant, is_constant_reused=is_constant_reused
            )
            if take_modulo
            else Script.parse_string(" ".join(["OP_FROMALTSTACK"] * (self.EXTENSION_DEGREE - 1)))
        )

        return out

    def square(
        self,
        take_modulo: bool,
        positive_modulo: bool = True,
        check_constant: bool | None = None,
        clean_constant: bool | None = None,
        is_constant_reused: bool | None = None,
        x: StackFiniteFieldElement | None = None,
        rolling_options: int = 1,
    ) -> Script:
        pass
