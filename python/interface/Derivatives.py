#!/usr/bin/env python3
"""
Wrapper class for easier printing of derivatives and functions.
"""

class Derivatives:
    """This is a nice container for the function and its derivatives, along with some
    useful helper functions.

    :param float f: function value in x
    :param float df_deta: ∂f/∂η
    :param float df_dtheta: ∂f/∂θ
    :param float d2f_deta2: ∂²f/∂η²
    :param float d2f_dtheta2: ∂²f/∂θ²
    :param float d2f_deta_dtheta: ∂²f/∂η∂θ
    """
    def __init__(self,*, f, df_deta, df_dtheta, d2f_deta2, d2f_dtheta2, d2f_deta_dtheta):
        self.f = f
        self.df_deta = df_deta
        self.df_dtheta = df_dtheta
        self.d2f_deta2 = d2f_deta2
        self.d2f_dtheta2 = d2f_dtheta2
        self.d2f_deta_dtheta = d2f_deta_dtheta

    def __str__(self):
        return f"""
 F       : {self.f},
 ∂f/∂η   : {self.df_deta},
 ∂f/∂θ   : {self.df_dtheta},
 ∂²f/∂η² : {self.d2f_deta2},
 ∂²f/∂θ² : {self.d2f_dtheta2},
 ∂²f/∂η∂θ: {self.d2f_deta_dtheta}"""

