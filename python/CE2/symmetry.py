import transpose
import reverse
def enforce_symmetry(field, symm_field=None):
    #logger.info("Enforcing symmetry.")
    Trans = transpose.TransposeOperator
    Rev = reverse.ReverseFirst

    if symm_field:
        symm = Rev(Trans(symm_field)).evaluate()
    else:
        # symmetry transform is transpose y1, y2, then reverse x
        symm = Rev(Trans(field)).evaluate()

    field['g'] = 0.5*(symm['g'] + field['g'])

    if symm_field:
        symm = Rev(Trans(field)).evaluate()
        symm_field['g'] = symm['g']
