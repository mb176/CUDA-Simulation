(* :Name: User`FixedTalbotNumericalLaplaceInversion` *)

(* :Title: Numerical Inversion of Laplace Transform with Multiple Precision Using the Complex Domain*)

(* :Authors: Joe Abate and Peter P. Valko*)

(* :Summary:
This package provides only one function: FT. The function calculates the value of the inverse of a Laplace transform at a specified time point. The Laplace transform should be provided as a function ready for multiple-precision evaluation in the complex plane. In other words, approximate numbers (with decimal point) or Mathematica functions starting with the letter 'N' are not allowed. The algorithm used is called Fixed Talbot algorithm.
*)

(* :Context: User`FixedTalbotNumericalLaplaceInversion` *)

(* :Package Version: 1.0 *)

(* :Copyright: Copyright 2003,  Joe Abate and Peter P. Valko *)

(* :History: Mathematica code originally written by Peter P. Valko, Sept 01, 2003. *)

(* :Keywords:
Laplace transform, Numerical inversion, Multiple-precision, Talbot method, Complex arithmetic 
*)

(* :Source:
Abate, J. and Valkó, P. P.:  Multi-Precision Laplace Transform Inversion, International Journal of Numerical Methods for Engineering, (2003), accepted for publication (IJNME A1032.C)
*)

(* :Warnings: None. *)

(* :Mathematica Version: 4.1 *)

(* :Limitations: 
The Laplace transform function should be given in a form, that permits multiple precision evaluation in the complex plane.
*)

(* :Discussion:
There are two main groups of methods for the numerical inversion of the Laplace transform. One of them uses values of the transform calculated on the real line. Such an algorithm (GWR) is available from  http://library.wolfram.com/infocenter/MathSource/4738/ . This submission falls into the other group of methods, requiring the evaluation of the transform in the complex plane. The best such method is due to Talbot and this realization gives an improved version of his algorithm. Deatils of the algorith are discussed in the paper Abate, J. and Valkó, P. P.:  Multi-Precision Laplace Transform Inversion, International Journal of Numerical Methods for Engineering, (2003), accepted for publication (IJNME A1032.C)
*)

BeginPackage["User`FixedTalbotNumericalLaplaceInversion`"]

FT::usage =
"FT[F, t, M] gives the inverse of the Laplace transform function named 'F' for a given time point 't'. The method involves the calculation of 'M' points on the modified Talbot path. The precision of internal calculations is set automatically.
\n
\n
FT[F, t] calculates the same but the number of points 'M' is set to the default 32. The precision of internal calculations is set automatically.
\n
\n
Important note: The Laplace transform should be defined as a function of one argument that can be complex number. It can involve anything from a simple Mathematica expression to a sophisticated Module. Since the Laplace transform will be evaluated with non-standard (multiple) precision, approximate numbers (with decimal point) or Mathematica functions starting with the letter 'N' are not allowed in the function definition.
\n
\n
Example usage:
\n
\n
fun[s_]=(1/s) Exp[-Sqrt[ (s^2 + (37/100) s + 1)/(s^2 + s + Pi)]]
\n
t0=100
\n
FT[fun,t0]"

Unprotect[FT];

Begin["User`FixedTalbotNumericalLaplaceTransformInversion`Private`"]
FT[F_, t_, M_:32]:=
    Module[{np, r, S, theta, sigma},
    np = Max[M, $MachinePrecision];
    r = SetPrecision[2M/(5t), np];
    S = r theta (Cot[theta] + I);
    sigma = theta + (theta Cot[theta] - 1)Cot[theta];
    (r/M)Plus @@ Append[Table[Re[Exp[t S](1 + I sigma)F[S]],
            {theta, Pi/M, (M - 1)Pi/M, Pi/M}],
                (1/2) Exp[r t] F[r]]
    ]
End[]  (* User`FixedTalbotNumericalLaplaceTransformInversion` *)
Protect[FT];
EndPackage[] (* User`FixedTalbotNumericalLaplaceInversion` *)
