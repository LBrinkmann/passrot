# passrot
Calculation of anisotropic scattering pattern.

Help:
passrot10.py -h

Example input:
passrot10.py -fa rerun/4B9O_pR0/?/waxs?_averageA.xvg -fb rerun/1TS8_pG/?/waxs?_averageA.xvg -env envelope/env_ft.dat -phi 360 -weight cos2 -polar linear -mtype linear -beam 12 -bulk 334 -legend "pR0" -t 10 -D "0.000021" -illustrate 0.1 1.5 3.0 4.5 6 -m -0.0284567 0.965887 0.257388 -o "output"
