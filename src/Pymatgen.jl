module Pymatgen
    using PythonCall


    export py_Struc


    const py_Struc = PythonCall.pynew()
    
    function __init__()
        PythonCall.pycopy!(
            py_Struc,
            pyimport("pymatgen.core").Structure
        )
    end
end