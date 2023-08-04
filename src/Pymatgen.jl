module Pymatgen
    using PythonCall


    export py_Struc, py_Incar
    export get_py_struc


    const py_Struc = PythonCall.pynew()
    const py_Incar = PythonCall.pynew()
    
    function __init__()
        PythonCall.pycopy!(
            py_Struc,
            pyimport("pymatgen.core").Structure
        )
        PythonCall.pycopy!(
            py_Incar,
            pyimport("pymatgen.io.vasp.inputs").Incar
        )
    end


    get_py_struc(filepath) = py_Struc.from_file(filepath)
end