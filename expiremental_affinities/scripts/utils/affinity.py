import math

def ktg(ki: float,temp_C: float,units:str = 'M') -> float:
    """
    Convert Ki/Kd in M to delta G in kcal/mol.
    Args:
        ki: Ki or Kd in Molar (M by default)
        temp_C: Temperature in Celsius
        units: Units of the input Ki/Kd. Can be 'M' for Molar, 'nM' for nanomolar, etc.
    Returns:
        Delta G in kcal/mol
    """

    R = 1.9872036e-3  # kcal/(mol*K)
    T = temp_C + 273.15  # convert Celsius to Kelvin
    # Convert Ki to Molar if necessary
    if units == 'M':
        ki_M = ki
    elif units == 'nM':
        ki_M = ki * 1e-9
        
    return R * T * math.log(ki_M)
