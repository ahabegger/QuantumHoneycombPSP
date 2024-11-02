from sympy import symbols, simplify, sympify

# Define variables
q1_1, q2_1, q2_2, q3_1, q3_2 = symbols('q1_1 q2_1 q2_2 q3_1 q3_2')

# Map for replacing variables in the expression
variable_map = {
    'q^1_1': 'q1_1',
    'q^2_1': 'q2_1',
    'q^2_2': 'q2_2',
    'q^3_1': 'q3_1',
    'q^3_2': 'q3_2'
}

# Preprocess the expression to replace variables and clean up syntax
for old_var, new_var in variable_map.items():
    expression = expression.replace(old_var, new_var)

# Remove newlines and extra spaces
expression = expression.replace('\n', ' ').replace('\t', ' ').strip()

# Since the expression is complex and may contain unsupported syntax for sympy,
# we'll need to further preprocess it to ensure it's valid.
# For this example, we'll assume the expression has been properly formatted.
# If necessary, you might need to perform additional replacements or parsing.

# Now, attempt to simplify the expression
try:
    # Convert the string expression to a sympy expression
    sympy_expr = sympify(expression)

    # Simplify the sympy expression
    simplified_expr = simplify(sympy_expr)

    # Print the simplified expression
    print("Simplified Energy Function:")
    print(simplified_expr)
except Exception as e:
    print("An error occurred while simplifying the expression:")
    print(e)
