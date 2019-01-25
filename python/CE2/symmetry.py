import transpose
def enforce_symmetry(field):
    #logger.info("Enforcing symmetry.")
    T = transpose.TransposeOperator
    
    trans = T(field).evaluate()

    field['g'] = 0.5*(trans['g'] + field['g'])
