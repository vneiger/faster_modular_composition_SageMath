####################
#  Some constants  #
####################

# Constant ETA which rules the choice of parameter m
omega = 2.37
omega_2 = 3.25
ETA = 1 / (1 + 2. * (omega-1) / (omega_2-2))

# Constants for failure and for certification
FAIL = "Fail"
CERT = "Cert"
NOCERT = "NoCert"

# Constant for debugging in examples.sage
DEBUG = False
