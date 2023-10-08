module Python
    using PythonCall


    export pyconvert
    export py_Struc, py_Incar, py_Outcar, py_Sga, py_np


    const py_Struc = PythonCall.pynew()
    const py_Incar = PythonCall.pynew()
    const py_Outcar = PythonCall.pynew()
    const py_Sga = PythonCall.pynew()
    const py_np = PythonCall.pynew()

    function __init__()
        PythonCall.pycopy!(
            py_Struc,
            pyimport("pymatgen.core").Structure
        )
        PythonCall.pycopy!(
            py_Incar,
            pyimport("pymatgen.io.vasp.inputs").Incar
        )
        PythonCall.pycopy!(
            py_Outcar,
            pyimport("pymatgen.io.vasp").Outcar
        )
        PythonCall.pycopy!(
            py_Sga,
            pyimport("pymatgen.symmetry.analyzer").SpacegroupAnalyzer
        )
        PythonCall.pycopy!(
            py_np,
            pyimport("numpy")
        )
    end
end