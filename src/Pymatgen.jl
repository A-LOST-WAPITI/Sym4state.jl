module Pymatgen
    using PythonCall


    export py_Struc, get_py_struc


    const py_Struc = PythonCall.pynew()
    
    function __init__()
        PythonCall.pycopy!(
            py_Struc,
            pyimport("pymatgen.core").Structure
        )
    end


    get_py_struc(filepath) = py_Struc.from_file(filepath)
end