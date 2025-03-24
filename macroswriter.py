def writeLatexMacro(macro_name:str, value, unit=None, error=None, digitsIfNoError:int = 4, filepath='Paper/macros.tex'):
    """
    Writes or overrides LaTeX macros in the specified file.
    Rounds the value and error to the same order of magnitude (two digits of the error) and writes them in scientific notation.
    """

    
    if error:
        value_order = int(f"{value:.1e}".split('e')[1])
        error_order = int(f"{error:.1e}".split('e')[1])
        error *= 10**-value_order
        value *= 10**-value_order
        valueStr = f"{value:.{value_order-error_order+1}f}"
        errorStr = f"{error:.{value_order -error_order+1}f}"

        macro_content = f"\\newcommand{{\\{macro_name}}}{{({valueStr} \\pm {errorStr}) \\times 10^{{{value_order}}}"
    else:
        value_order = int(f"{value:.1e}".split('e')[1])
        value *= 10**-value_order
        valueStr = f"{value:.{digitsIfNoError-1}f}"

        macro_content = f"\\newcommand{{\\{macro_name}}}{{{valueStr} \\times 10^{{{value_order}}}"

    # Add unit if given
    if unit:
        macro_content += f"\\,{unit}"
    macro_content += "}\n"

    # Read the existing file content
    try:
        with open(filepath, 'r') as file:
            lines = file.readlines()
    except FileNotFoundError:
        lines = []

    # Remove existing macro with the same name
    lines = [line for line in lines if not line.strip().startswith(f"\\newcommand{{\\{macro_name}}}")] 

    # Append the new macro
    lines.append(macro_content)

    # Write back to the file
    with open(filepath, 'w') as file:
        file.writelines(lines)





