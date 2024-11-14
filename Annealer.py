import pickle
import time
from dimod import ConstrainedQuadraticModel, Binary

def annealer(bqm):
    cqm = ConstrainedQuadraticModel().from_bqm(bqm)

    linear_variables = []
    ancillary_variables = []
    variables = []

    print(bqm.to_polystring())
    print(bqm.variables)

    for var in bqm.variables:
        variables.append(var)
        if '*' in var:
            ancillary_variables.append(var)
        else:
            linear_variables.append(var)

    binary_variables = {var: Binary(var) for var in variables}

    for ancillary in ancillary_variables:
        stripped_bit = ancillary.replace(' ', '')
        base_bits = stripped_bit.split('*')
        a = binary_variables[ancillary]
        for bit in base_bits:
            b = binary_variables[bit]
            cqm.add_constraint(a - b <= 0, label=f'{ancillary} >= {bit}')

        sum_base_bits = sum(binary_variables[base_bit] for base_bit in base_bits)
        length_constraint = -1 * (len(base_bits) - 1)
        cqm.add_constraint(a - sum_base_bits >= length_constraint, label=f'{ancillary} >= {base_bits}')

    # Print all constraints
    for label, constraint in cqm.constraints.items():
        print(f"Constraint '{label}': {constraint.lhs} {constraint.sense} {constraint.rhs}")

    api_token = '' # Insert your API token here

    from dwave.system import LeapHybridCQMSampler
    sampler = LeapHybridCQMSampler(token=api_token)
    sampleset = sampler.sample_cqm(cqm)

    filtered_samples = sampleset.filter(lambda row: row.is_feasible)
    filtered_samples = filtered_samples.aggregate()

    print(filtered_samples.record.sample)
    print(filtered_samples.record.num_occurrences)
    print(filtered_samples.record.energy)

    # Pickle the filtered samples
    with open(f'Results/Samples_{int(time.time())}.pkl', 'wb') as f:
        pickle.dump(filtered_samples, f)

    return filtered_samples
