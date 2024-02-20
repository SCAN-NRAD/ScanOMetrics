"""
Miscellaneous util functions
"""

import random
import string

def _generate_random_code(length=10, seed=None):
    """Generate random alphanumerical code, starting with a letter for proper sorting in file browsers.

    :param length: specifies length of random string, defaults to 10
    :type length: int
    :param seed: seed for random string generator, defaults to None
    :type seed: int, optional
    :return: random alphanumerical string with given length
    :rtype: string
    """

    random.seed(seed)
    return ''.join(random.choices(string.ascii_uppercase)+random.choices(string.ascii_uppercase+string.digits,k=length-1))

def _generate_safeReading_random_code(length=10, seed=None):
    """Generate random alphanumerical code from a set without commonly misread characters
    (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3541865), starting with a letter for proper sorting in file browsers.

    :param length: specifies length of random string, defaults to 10
    :type length: int
    :param seed: seed for random string generator, defaults to None
    :type seed: int, optional
    :return: random alphanumerical string with given length
    :rtype: string
    """

    safe_alpha_set = 'ACGHJKLMNPQRUVWXY'
    safe_alphanumerical_set = 'ACGHJKLMNPQRUVWXY469'
    random.seed(seed)
    return ''.join(random.choices(safe_alpha_set)+random.choices(safe_alphanumerical_set, k=length-1))
