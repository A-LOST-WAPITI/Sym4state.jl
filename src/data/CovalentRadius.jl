# Data from https://doi.org/10.1039/B801115J

const RADIUS_DICT = Dict(
    "H" => Dict(
        "name" => ["H"],
        "radii" => [0.31]
    ),
    "He" => Dict(
        "name" => ["He"],
        "radii" => [0.28]
    ),
    "Li" => Dict(
        "name" => ["Li"],
        "radii" => [1.28]
    ),
    "Be" => Dict(
        "name" => ["Be"],
        "radii" => [0.96]
    ),
    "B" => Dict(
        "name" => ["B"],
        "radii" => [0.84]
    ),
    "C" => Dict(
        "name" => ["Csp3", "Csp2", "Csp"],
        "radii" => [0.76, 0.73, 0.69]
    ),
    "N" => Dict(
        "name" => ["N"],
        "radii" => [0.71]
    ),
    "O" => Dict(
        "name" => ["O"],
        "radii" => [0.66]
    ),
    "F" => Dict(
        "name" => ["F"],
        "radii" => [0.57]
    ),
    "Ne" => Dict(
        "name" => ["Ne"],
        "radii" => [0.58]
    ),
    "Na" => Dict(
        "name" => ["Na"],
        "radii" => [1.66]
    ),
    "Mg" => Dict(
        "name" => ["Mg"],
        "radii" => [1.41]
    ),
    "Al" => Dict(
        "name" => ["Al"],
        "radii" => [1.21]
    ),
    "Si" => Dict(
        "name" => ["Si"],
        "radii" => [1.11]
    ),
    "P" => Dict(
        "name" => ["P"],
        "radii" => [1.07]
    ),
    "S" => Dict(
        "name" => ["S"],
        "radii" => [1.05]
    ),
    "Cl" => Dict(
        "name" => ["Cl"],
        "radii" => [1.02]
    ),
    "Ar" => Dict(
        "name" => ["Ar"],
        "radii" => [1.06]
    ),
    "K" => Dict(
        "name" => ["K"],
        "radii" => [2.03]
    ),
    "Ca" => Dict(
        "name" => ["Ca"],
        "radii" => [1.76]
    ),
    "Sc" => Dict(
        "name" => ["Sc"],
        "radii" => [1.7]
    ),
    "Ti" => Dict(
        "name" => ["Ti"],
        "radii" => [1.6]
    ),
    "V" => Dict(
        "name" => ["V"],
        "radii" => [1.53]
    ),
    "Cr" => Dict(
        "name" => ["Cr"],
        "radii" => [1.39]
    ),
    "Mn" => Dict(
        "name" => ["Mn", "Mn"],
        "radii" => SubString{String}["l.s.", "h.s."]
    ),
    "Fe" => Dict(
        "name" => ["Fe", "Fe"],
        "radii" => SubString{String}["l.s.", "h.s."]
    ),
    "Co" => Dict(
        "name" => ["Co", "Co"],
        "radii" => SubString{String}["l.s.", "h.s."]
    ),
    "Ni" => Dict(
        "name" => ["Ni"],
        "radii" => [1.24]
    ),
    "Cu" => Dict(
        "name" => ["Cu"],
        "radii" => [1.32]
    ),
    "Zn" => Dict(
        "name" => ["Zn"],
        "radii" => [1.22]
    ),
    "Ga" => Dict(
        "name" => ["Ga"],
        "radii" => [1.22]
    ),
    "Ge" => Dict(
        "name" => ["Ge"],
        "radii" => [1.2]
    ),
    "As" => Dict(
        "name" => ["As"],
        "radii" => [1.19]
    ),
    "Se" => Dict(
        "name" => ["Se"],
        "radii" => [1.2]
    ),
    "Br" => Dict(
        "name" => ["Br"],
        "radii" => [1.2]
    ),
    "Kr" => Dict(
        "name" => ["Kr"],
        "radii" => [1.16]
    ),
    "Rb" => Dict(
        "name" => ["Rb"],
        "radii" => [2.2]
    ),
    "Sr" => Dict(
        "name" => ["Sr"],
        "radii" => [1.95]
    ),
    "Y" => Dict(
        "name" => ["Y"],
        "radii" => [1.9]
    ),
    "Zr" => Dict(
        "name" => ["Zr"],
        "radii" => [1.75]
    ),
    "Nb" => Dict(
        "name" => ["Nb"],
        "radii" => [1.64]
    ),
    "Mo" => Dict(
        "name" => ["Mo"],
        "radii" => [1.54]
    ),
    "Tc" => Dict(
        "name" => ["Tc"],
        "radii" => [1.47]
    ),
    "Ru" => Dict(
        "name" => ["Ru"],
        "radii" => [1.46]
    ),
    "Rh" => Dict(
        "name" => ["Rh"],
        "radii" => [1.42]
    ),
    "Pd" => Dict(
        "name" => ["Pd"],
        "radii" => [1.39]
    ),
    "Ag" => Dict(
        "name" => ["Ag"],
        "radii" => [1.45]
    ),
    "Cd" => Dict(
        "name" => ["Cd"],
        "radii" => [1.44]
    ),
    "In" => Dict(
        "name" => ["In"],
        "radii" => [1.42]
    ),
    "Sn" => Dict(
        "name" => ["Sn"],
        "radii" => [1.39]
    ),
    "Sb" => Dict(
        "name" => ["Sb"],
        "radii" => [1.39]
    ),
    "Te" => Dict(
        "name" => ["Te"],
        "radii" => [1.38]
    ),
    "I" => Dict(
        "name" => ["I"],
        "radii" => [1.39]
    ),
    "Xe" => Dict(
        "name" => ["Xe"],
        "radii" => [1.4]
    ),
    "Cs" => Dict(
        "name" => ["Cs"],
        "radii" => [2.44]
    ),
    "Ba" => Dict(
        "name" => ["Ba"],
        "radii" => [2.15]
    ),
    "La" => Dict(
        "name" => ["La"],
        "radii" => [2.07]
    ),
    "Ce" => Dict(
        "name" => ["Ce"],
        "radii" => [2.04]
    ),
    "Pr" => Dict(
        "name" => ["Pr"],
        "radii" => [2.03]
    ),
    "Nd" => Dict(
        "name" => ["Nd"],
        "radii" => [2.01]
    ),
    "Pm" => Dict(
        "name" => ["Pm"],
        "radii" => [1.99]
    ),
    "Sm" => Dict(
        "name" => ["Sm"],
        "radii" => [1.98]
    ),
    "Eu" => Dict(
        "name" => ["Eu"],
        "radii" => [1.98]
    ),
    "Gd" => Dict(
        "name" => ["Gd"],
        "radii" => [1.96]
    ),
    "Tb" => Dict(
        "name" => ["Tb"],
        "radii" => [1.94]
    ),
    "Dy" => Dict(
        "name" => ["Dy"],
        "radii" => [1.92]
    ),
    "Ho" => Dict(
        "name" => ["Ho"],
        "radii" => [1.92]
    ),
    "Er" => Dict(
        "name" => ["Er"],
        "radii" => [1.89]
    ),
    "Tm" => Dict(
        "name" => ["Tm"],
        "radii" => [1.9]
    ),
    "Yb" => Dict(
        "name" => ["Yb"],
        "radii" => [1.87]
    ),
    "Lu" => Dict(
        "name" => ["Lu"],
        "radii" => [1.87]
    ),
    "Hf" => Dict(
        "name" => ["Hf"],
        "radii" => [1.75]
    ),
    "Ta" => Dict(
        "name" => ["Ta"],
        "radii" => [1.7]
    ),
    "W" => Dict(
        "name" => ["W"],
        "radii" => [1.62]
    ),
    "Re" => Dict(
        "name" => ["Re"],
        "radii" => [1.51]
    ),
    "Os" => Dict(
        "name" => ["Os"],
        "radii" => [1.44]
    ),
    "Ir" => Dict(
        "name" => ["Ir"],
        "radii" => [1.41]
    ),
    "Pt" => Dict(
        "name" => ["Pt"],
        "radii" => [1.36]
    ),
    "Au" => Dict(
        "name" => ["Au"],
        "radii" => [1.36]
    ),
    "Hg" => Dict(
        "name" => ["Hg"],
        "radii" => [1.32]
    ),
    "Tl" => Dict(
        "name" => ["Tl"],
        "radii" => [1.45]
    ),
    "Pb" => Dict(
        "name" => ["Pb"],
        "radii" => [1.46]
    ),
    "Bi" => Dict(
        "name" => ["Bi"],
        "radii" => [1.48]
    ),
    "Po" => Dict(
        "name" => ["Po"],
        "radii" => [1.4]
    ),
    "At" => Dict(
        "name" => ["At"],
        "radii" => [1.5]
    ),
    "Rn" => Dict(
        "name" => ["Rn"],
        "radii" => [1.5]
    ),
    "Fr" => Dict(
        "name" => ["Fr"],
        "radii" => [2.6]
    ),
    "Ra" => Dict(
        "name" => ["Ra"],
        "radii" => [2.21]
    ),
    "Ac" => Dict(
        "name" => ["Ac"],
        "radii" => [2.15]
    ),
    "Th" => Dict(
        "name" => ["Th"],
        "radii" => [2.06]
    ),
    "Pa" => Dict(
        "name" => ["Pa"],
        "radii" => [2.0]
    ),
    "U" => Dict(
        "name" => ["U"],
        "radii" => [1.96]
    ),
    "Np" => Dict(
        "name" => ["Np"],
        "radii" => [1.9]
    ),
    "Pu" => Dict(
        "name" => ["Pu"],
        "radii" => [1.87]
    ),
    "Am" => Dict(
        "name" => ["Am"],
        "radii" => [1.8]
    ),
    "Cm" => Dict(
        "name" => ["Cm"],
        "radii" => [1.69]
    ),
)