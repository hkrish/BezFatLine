
var EPSILON = 10e-12;
var TOLERANCE = 10e-6;
var MAX_ITERATE = 20;

// :)
var BinomialCoeffs = [
                        [1],
                        [1, 1],
                        [1, 2, 1],
                        [1, 3, 3, 1],
                        [1, 4, 6, 4, 1],
                        [1, 5, 10, 10, 5, 1],
                        [1, 6, 15, 20, 15, 6, 1],
                        [1, 7, 21, 35, 35, 21, 7, 1],
                        [1, 8, 28, 56, 70, 56, 28, 8, 1],
                        [1, 9, 36, 84, 126, 126, 84, 36, 9, 1]
                    ];

function evaluateHorner( v, t ) {
    var i, n = v.length-1, h = v[n];
    for( i = n-1; i >= 0; i-- )
        h = t * h + v[i];
    return h;
}

// v[0] is the coefficient for term t^0, v[1] - t^1 .. v[n] - t^n ...
function polyLongDivide( v, d ){
    var n = v.length-1, h = v[n], dh, i;
    var v1 = [];
    for( i = n-1; i >= 0; i-- ){
        dh = d * h;
        v1.unshift( dh );
        h= dh + v[i];
    }
    return v1;
}

// v[0] is the coefficient for term t^0, v[1] - t^1 .. v[n] - t^n ...
function polyLongDivideSelf( v, d ){
    var n = v.length-1, h = v[n], dh, i, t;
    for( i = n-1; i >= 0; i-- ){
        dh = d * h;
        t = v[i];
        v[i] = dh / d;
        h= dh + t;
    }
    v.pop();
}

// returns the power basis for for the bezier curvre of degree n
// vp[0] is the coefficient for term t^0, vp[1] - t^1 .. vp[n] - t^n ...
function bernsteinToBasis( vb ){
    var n = vb.length-1, vp = [], vd = [], i, j, l;
    var bcoeffs = BinomialCoeffs[ n ];
    vp.push( vb[0] );
    for( i = 0; i <= n; i++ )
        vd[i] = vb[i];
    for( i = 1; i <= n; i++ ) {
        for( j = 0, l = n-i; j <= l; j++ ) {
            vd[j] = vd[j+1] - vd[j];
        }
        vp.push( vd[0] * bcoeffs[i] );
    }
    return vp;
}

// For Cubic bezier curves
//
// References:
// [SP86] - Thomas W. Sederberg and Scott R. Parry.
//      A comparison of curve-curve intersection algorithms.
//      Computer-Aided Design, 18:58–63, 1986.
function bernsteinToImplicit( v ){
    var x1 = v[0], y1 = v[1], x2 = v[2], y2 = v[3];
    var x3 = v[4], y3 = v[5], x4 = v[6], y4 = v[7];
    // Squares and cubes of coefficients
    var x12 = x1*x1, y12 = y1*y1, x22 = x2*x2, y22 = y2*y2;
    var x32 = x3*x3, y32 = y3*y3, x42 = x4*x4, y42 = y4*y4;
    var x13= x12*x1, y13= y12*y1, x23= x22*x2, y23= y22*y2;
    var x33= x32*x3, y33= y32*y3, x43= x42*x4, y43= y42*y4;

    /*
    * According to [SP86], the implicit form of a cubic bezier curve
    * from its parametric form is
    *
    *           |  l32      l31       l30  |
    * f(x, y) = |  l31   l30 + l21    l20  |
    *           |  l30      l20       l10  |
    *
    *  where,
    *                         | x   y   1 |
    *      lij = (3 i)*(3 j)* | xi  yi  1 |
    *                         | xj  yj  1 |
    *
    *        * Note: (3 i) is the binomial coeff. 3! / ( i! * (3-i)! )
    *
    *  let P=l32, Q=l31, R=l30, S=l21, T=l20 and U=l10
    * f(x, y) = P*R*U + P*S*U - P*(T^2) - (Q^2)*U + 2*Q*R*T - (R^3) - (R^2)*S
    *
    * Implicit form is returned as an array
    *      [ A, B, C, D, E, F, G, H, I, J ]
    * for the implicit form
    *      F(x, y) = A y^3 + B y^2 + C y + D x^3 + x^2 ( E + F y ) + x ( G + H y^2 + I y ) + J
    *
    * Expansion of these coefficients is as follows (thanks to Mathematica!)
    */

    var A = -x13+9*x12*x2-27*x1*x22+27*x23-9*x12*x3+54*x1*x2*x3-
        81*x22*x3-27*x1*x32+81*x2*x32-27*x33+3*x12*x4-
        18*x1*x2*x4+27*x22*x4+18*x1*x3*x4-54*x2*x3*x4+
        27*x32*x4-3*x1*x42+9*x2*x42-9*x3*x42+x43;

    var B = -27*x23*y1+27*x1*x2*x3*y1+81*x22*x3*y1-54*x1*x32*y1-
        81*x2*x32*y1+54*x33*y1-3*x12*x4*y1-9*x1*x2*x4*y1-
        27*x22*x4*y1+63*x1*x3*x4*y1+27*x2*x3*x4*y1-54*x32*x4*y1-
        21*x1*x42*y1+9*x2*x42*y1+18*x3*x42*y1-3*x43*y1+
        27*x1*x22*y2-18*x12*x3*y2-81*x1*x2*x3*y2+162*x1*x32*y2-
        81*x2*x32*y2+27*x12*x4*y2-27*x1*x2*x4*y2+54*x22*x4*y2-
        153*x1*x3*x4*y2+81*x2*x3*x4*y2+54*x1*x42*y2-
        54*x2*x42*y2+9*x3*x42*y2-9*x12*x2*y3+54*x12*x3*y3-
        81*x1*x2*x3*y3+81*x22*x3*y3-54*x1*x32*y3-54*x12*x4*y3+
        153*x1*x2*x4*y3-162*x22*x4*y3+27*x1*x3*x4*y3+
        81*x2*x3*x4*y3-27*x32*x4*y3-27*x1*x42*y3+18*x2*x42*y3+
        3*x13*y4-18*x12*x2*y4+54*x1*x22*y4-54*x23*y4-
        9*x12*x3*y4-27*x1*x2*x3*y4+81*x22*x3*y4+27*x1*x32*y4-
        81*x2*x32*y4+27*x33*y4+21*x12*x4*y4-63*x1*x2*x4*y4+
        54*x22*x4*y4+9*x1*x3*x4*y4-27*x2*x3*x4*y4+3*x1*x42*y4;


    var C = -27*x33*y12+27*x2*x3*x4*y12+27*x32*x4*y12-
        3*x1*x42*y12-18*x2*x42*y12-9*x3*x42*y12+3*x43*y12+
        81*x2*x32*y1*y2-54*x22*x4*y1*y2-9*x1*x3*x4*y1*y2-
        81*x2*x3*x4*y1*y2+27*x1*x42*y1*y2+54*x2*x42*y1*y2-
        18*x3*x42*y1*y2-81*x1*x32*y22+54*x1*x2*x4*y22+
        81*x1*x3*x4*y22-81*x1*x42*y22+27*x2*x42*y22-
        81*x22*x3*y1*y3+54*x1*x32*y1*y3+9*x1*x2*x4*y1*y3+
        162*x22*x4*y1*y3-108*x1*x3*x4*y1*y3-81*x2*x3*x4*y1*y3+
        54*x32*x4*y1*y3+27*x1*x42*y1*y3-36*x2*x42*y1*y3+
        81*x1*x2*x3*y2*y3-27*x12*x4*y2*y3-162*x1*x2*x4*y2*y3+
        162*x1*x3*x4*y2*y3-81*x2*x3*x4*y2*y3+27*x1*x42*y2*y3-
        27*x12*x3*y32+81*x12*x4*y32-81*x1*x2*x4*y32+
        81*x22*x4*y32-54*x1*x3*x4*y32+54*x23*y1*y4-
        54*x1*x2*x3*y1*y4-81*x22*x3*y1*y4+54*x1*x32*y1*y4+
        81*x2*x32*y1*y4-54*x33*y1*y4+6*x12*x4*y1*y4+
        9*x1*x2*x4*y1*y4-54*x22*x4*y1*y4-9*x1*x3*x4*y1*y4+
        54*x2*x3*x4*y1*y4-6*x1*x42*y1*y4-54*x1*x22*y2*y4+
        36*x12*x3*y2*y4+81*x1*x2*x3*y2*y4-162*x1*x32*y2*y4+
        81*x2*x32*y2*y4-27*x12*x4*y2*y4+108*x1*x2*x4*y2*y4-
        54*x22*x4*y2*y4-9*x1*x3*x4*y2*y4+18*x12*x2*y3*y4-
        54*x12*x3*y3*y4+81*x1*x2*x3*y3*y4-81*x22*x3*y3*y4+
        54*x1*x32*y3*y4-27*x12*x4*y3*y4+9*x1*x2*x4*y3*y4-
        3*x13*y42+9*x12*x2*y42-27*x1*x22*y42+27*x23*y42+
        18*x12*x3*y42-27*x1*x2*x3*y42+3*x12*x4*y42;

    var D = y13-9*y12*y2+27*y1*y22-27*y23+9*y12*y3-
        54*y1*y2*y3+81*y22*y3+27*y1*y32-81*y2*y32+27*y33-
        3*y12*y4+18*y1*y2*y4-27*y22*y4-18*y1*y3*y4+
        54*y2*y3*y4-27*y32*y4+3*y1*y42-9*y2*y42+9*y3*y42-
        y43;

    var E = -3*x4*y13+9*x3*y12*y2+18*x4*y12*y2-
        27*x2*y1*y22-54*x4*y1*y22+27*x1*y23+54*x4*y23+
        18*x2*y12*y3-54*x3*y12*y3+9*x4*y12*y3-27*x1*y1*y2*y3+
        81*x2*y1*y2*y3+81*x3*y1*y2*y3+27*x4*y1*y2*y3-
        81*x1*y22*y3-81*x3*y22*y3-81*x4*y22*y3+54*x1*y1*y32-
        162*x2*y1*y32+54*x3*y1*y32-27*x4*y1*y32+81*x1*y2*y32+
        81*x2*y2*y32+81*x4*y2*y32-54*x1*y33-27*x4*y33+
        3*x1*y12*y4-27*x2*y12*y4+54*x3*y12*y4-21*x4*y12*y4+
        9*x1*y1*y2*y4+27*x2*y1*y2*y4-153*x3*y1*y2*y4+
        63*x4*y1*y2*y4+27*x1*y22*y4-54*x2*y22*y4+162*x3*y22*y4-
        54*x4*y22*y4-63*x1*y1*y3*y4+153*x2*y1*y3*y4-
        27*x3*y1*y3*y4-9*x4*y1*y3*y4-27*x1*y2*y3*y4-
        81*x2*y2*y3*y4-81*x3*y2*y3*y4+27*x4*y2*y3*y4+
        54*x1*y32*y4+27*x3*y32*y4+21*x1*y1*y42-54*x2*y1*y42+
        27*x3*y1*y42-3*x4*y1*y42-9*x1*y2*y42+54*x2*y2*y42-
        18*x3*y2*y42-18*x1*y3*y42-9*x2*y3*y42+3*x1*y43;

    var F = -3*x1*y12+9*x2*y12-9*x3*y12+3*x4*y12+18*x1*y1*y2-
        54*x2*y1*y2+54*x3*y1*y2-18*x4*y1*y2-27*x1*y22+
        81*x2*y22-81*x3*y22+27*x4*y22-18*x1*y1*y3+
        54*x2*y1*y3-54*x3*y1*y3+18*x4*y1*y3+54*x1*y2*y3-
        162*x2*y2*y3+162*x3*y2*y3-54*x4*y2*y3-27*x1*y32+
        81*x2*y32-81*x3*y32+27*x4*y32+6*x1*y1*y4-
        18*x2*y1*y4+18*x3*y1*y4-6*x4*y1*y4-18*x1*y2*y4+
        54*x2*y2*y4-54*x3*y2*y4+18*x4*y2*y4+18*x1*y3*y4-
        54*x2*y3*y4+54*x3*y3*y4-18*x4*y3*y4-3*x1*y42+
        9*x2*y42-9*x3*y42+3*x4*y42;

    var G = 3*x42*y13-18*x3*x4*y12*y2-9*x42*y12*y2+
        54*x2*x4*y1*y22+27*x42*y1*y22-54*x1*x4*y23-
        27*x42*y23+27*x32*y12*y3-36*x2*x4*y12*y3+
        54*x3*x4*y12*y3-18*x42*y12*y3-81*x2*x3*y1*y2*y3+
        54*x1*x4*y1*y2*y3-81*x2*x4*y1*y2*y3-81*x3*x4*y1*y2*y3+
        27*x42*y1*y2*y3+81*x1*x3*y22*y3+81*x1*x4*y22*y3+
        81*x3*x4*y22*y3+81*x22*y1*y32-54*x1*x3*y1*y32-
        54*x1*x4*y1*y32+162*x2*x4*y1*y32-54*x3*x4*y1*y32-
        81*x1*x2*y2*y32-81*x1*x4*y2*y32-81*x2*x4*y2*y32+
        27*x12*y33+54*x1*x4*y33+27*x2*x3*y12*y4-
        81*x32*y12*y4-6*x1*x4*y12*y4+27*x2*x4*y12*y4+
        27*x3*x4*y12*y4-3*x42*y12*y4-54*x22*y1*y2*y4-
        9*x1*x3*y1*y2*y4+162*x2*x3*y1*y2*y4+81*x32*y1*y2*y4-
        9*x1*x4*y1*y2*y4-108*x2*x4*y1*y2*y4-9*x3*x4*y1*y2*y4+
        54*x1*x2*y22*y4-162*x1*x3*y22*y4-81*x32*y22*y4+
        54*x1*x4*y22*y4+54*x2*x4*y22*y4+9*x1*x2*y1*y3*y4-
        81*x22*y1*y3*y4+108*x1*x3*y1*y3*y4-162*x2*x3*y1*y3*y4+
        54*x32*y1*y3*y4+9*x1*x4*y1*y3*y4+9*x2*x4*y1*y3*y4-
        27*x12*y2*y3*y4+81*x1*x2*y2*y3*y4+81*x1*x3*y2*y3*y4+
        81*x2*x3*y2*y3*y4-54*x1*x4*y2*y3*y4-27*x12*y32*y4-
        54*x1*x3*y32*y4+3*x12*y1*y42-27*x1*x2*y1*y42+
        81*x22*y1*y42-27*x1*x3*y1*y42-27*x2*x3*y1*y42+
        6*x1*x4*y1*y42+18*x12*y2*y42-54*x1*x2*y2*y42-
        27*x22*y2*y42+36*x1*x3*y2*y42+9*x12*y3*y42+
        18*x1*x2*y3*y42-3*x12*y43;

    var H = 3*x12*y1-18*x1*x2*y1+27*x22*y1+18*x1*x3*y1-
        54*x2*x3*y1+27*x32*y1-6*x1*x4*y1+18*x2*x4*y1-
        18*x3*x4*y1+3*x42*y1-9*x12*y2+54*x1*x2*y2-
        81*x22*y2-54*x1*x3*y2+162*x2*x3*y2-81*x32*y2+
        18*x1*x4*y2-54*x2*x4*y2+54*x3*x4*y2-9*x42*y2+
        9*x12*y3-54*x1*x2*y3+81*x22*y3+54*x1*x3*y3-
        162*x2*x3*y3+81*x32*y3-18*x1*x4*y3+54*x2*x4*y3-
        54*x3*x4*y3+9*x42*y3-3*x12*y4+18*x1*x2*y4-
        27*x22*y4-18*x1*x3*y4+54*x2*x3*y4-27*x32*y4+
        6*x1*x4*y4-18*x2*x4*y4+18*x3*x4*y4-3*x42*y4;

    var I = -27*x2*x3*y12+54*x32*y12+6*x1*x4*y12+9*x2*x4*y12-
        63*x3*x4*y12+21*x42*y12+54*x22*y1*y2+9*x1*x3*y1*y2-
        81*x2*x3*y1*y2-81*x32*y1*y2-45*x1*x4*y1*y2+
        81*x2*x4*y1*y2+126*x3*x4*y1*y2-63*x42*y1*y2-
        54*x1*x2*y22+81*x1*x3*y22+81*x32*y22+27*x1*x4*y22-
        108*x2*x4*y22-81*x3*x4*y22+54*x42*y22-
        9*x1*x2*y1*y3-81*x22*y1*y3+243*x2*x3*y1*y3-
        108*x32*y1*y3+45*x1*x4*y1*y3-180*x2*x4*y1*y3+
        81*x3*x4*y1*y3+9*x42*y1*y3+27*x12*y2*y3+
        81*x1*x2*y2*y3-243*x1*x3*y2*y3+243*x2*x4*y2*y3-
        81*x3*x4*y2*y3-27*x42*y2*y3-54*x12*y32+
        81*x1*x2*y32-81*x22*y32+108*x1*x3*y32-
        27*x1*x4*y32-81*x2*x4*y32+54*x3*x4*y32-6*x12*y1*y4+
        45*x1*x2*y1*y4-27*x22*y1*y4-45*x1*x3*y1*y4+
        27*x32*y1*y4+45*x2*x4*y1*y4-45*x3*x4*y1*y4+
        6*x42*y1*y4-9*x12*y2*y4-81*x1*x2*y2*y4+
        108*x22*y2*y4+180*x1*x3*y2*y4-243*x2*x3*y2*y4+
        81*x32*y2*y4-45*x1*x4*y2*y4+9*x3*x4*y2*y4+
        63*x12*y3*y4-126*x1*x2*y3*y4+81*x22*y3*y4-
        81*x1*x3*y3*y4+81*x2*x3*y3*y4-54*x32*y3*y4+
        45*x1*x4*y3*y4-9*x2*x4*y3*y4-21*x12*y42+
        63*x1*x2*y42-54*x22*y42-9*x1*x3*y42+27*x2*x3*y42-
        6*x1*x4*y42;

    var J = -x43*y13+9*x3*x42*y12*y2-27*x2*x42*y1*y22+
        27*x1*x42*y23-27*x32*x4*y12*y3+18*x2*x42*y12*y3+
        81*x2*x3*x4*y1*y2*y3-27*x1*x42*y1*y2*y3-81*x1*x3*x4*y22*y3-
        81*x22*x4*y1*y32+54*x1*x3*x4*y1*y32+81*x1*x2*x4*y2*y32-
        27*x12*x4*y33+27*x33*y12*y4-27*x2*x3*x4*y12*y4+
        3*x1*x42*y12*y4-81*x2*x32*y1*y2*y4+54*x22*x4*y1*y2*y4+
        9*x1*x3*x4*y1*y2*y4+81*x1*x32*y22*y4-54*x1*x2*x4*y22*y4+
        81*x22*x3*y1*y3*y4-54*x1*x32*y1*y3*y4-9*x1*x2*x4*y1*y3*y4-
        81*x1*x2*x3*y2*y3*y4+27*x12*x4*y2*y3*y4+27*x12*x3*y32*y4-
        27*x23*y1*y42+27*x1*x2*x3*y1*y42-3*x12*x4*y1*y42+
        27*x1*x22*y2*y42-18*x12*x3*y2*y42-9*x12*x2*y3*y42+
        x13*y43;

    var min = Math.min( A, B, C, D, E, F, G, H, I, J );
    return [ A/min, B/min, C/min, D/min, E/min, F/min, G/min, H/min, I/min, J/min ];
}

function getIntersectEquation( im, v ){
    var x1 = v[0], y1 = v[1], x2 = v[2], y2 = v[3];
    var x3 = v[4], y3 = v[5], x4 = v[6], y4 = v[7];
    // Squares and cubes of coefficients
    var x12 = x1*x1, y12 = y1*y1, x22 = x2*x2, y22 = y2*y2;
    var x32 = x3*x3, y32 = y3*y3, x42 = x4*x4, y42 = y4*y4;
    var x13= x12*x1, y13= y12*y1, x23= x22*x2, y23= y22*y2;
    var x33= x32*x3, y33= y32*y3, x43= x42*x4, y43= y42*y4;
    // Coefficients of the implicit equation
    var A = im[0], B = im[1], C = im[2], D = im[3];
    var E = im[4], F = im[5], G = im[6], H = im[7];
    var I = im[8], J = im[9];

    /*
     * One curve is of the implicit form
     *     f(x, y)
     * and another is given by it's parametric euation
     *     x = (1-t)^3 x1 + 3 (1-t)^2 t x2 + 3 (1-t) t^2 x3 + t^3 x4
     *     y = (1-t)^3 y1 + 3 (1-t)^2 t y2 + 3 (1-t) t^2 y3 + t^3 y4
     *
     * Substituting for x and y in the implicit equation f(x, y), we would
     * get a 9 degree polynomial, roots can be calculated by solving
     * this polynomial for real roots in t = [0..1]
     *
     * t0 is the coefficient for term t^0, t1 - t^1 .. t9 - t^9
     */

    var t0 = J+G*x1+E*x12+D*x13+C*y1+I*x1*y1+F*x12*y1+
        B*y12+H*x1*y12+A*y13;

    var t1 = G*(-3*x1+3*x2)+E*(-6*x12+6*x1*x2)+
        D*(-9*x13+9*x12*x2)+C*(-3*y1+3*y2)+
        I*(-6*x1*y1+3*x2*y1+3*x1*y2)+
        F*(-9*x12*y1+6*x1*x2*y1+3*x12*y2)+
        B*(-6*y12+6*y1*y2)+
        H*(-9*x1*y12+3*x2*y12+6*x1*y1*y2)+
        A*(-9*y13+9*y12*y2);

    var t2 = G*(3*x1-6*x2+3*x3)+
        E*(15*x12-30*x1*x2+9*x22+6*x1*x3)+
        D*(36*x13-72*x12*x2+27*x1*x22+9*x12*x3)+
        C*(3*y1-6*y2+3*y3)+
        I*(15*x1*y1-15*x2*y1+3*x3*y1-15*x1*y2+9*x2*y2+
        3*x1*y3)+
        F*(36*x12*y1-48*x1*x2*y1+9*x22*y1+6*x1*x3*y1-
        24*x12*y2+18*x1*x2*y2+3*x12*y3)+
        B*(15*y12-30*y1*y2+9*y22+6*y1*y3)+
        H*(36*x1*y12-24*x2*y12+3*x3*y12-48*x1*y1*y2+
        18*x2*y1*y2+9*x1*y22+6*x1*y1*y3)+
        A*(36*y13-72*y12*y2+27*y1*y22+9*y12*y3);

    var t3 = G*(-x1+3*x2-3*x3+x4)+
        E*(-20*x12+60*x1*x2-36*x22-24*x1*x3+18*x2*x3+
        2*x1*x4)+
        D*(-84*x13+252*x12*x2-189*x1*x22+27*x23-63*x12*x3+
        54*x1*x2*x3+3*x12*x4)+C*(-y1+3*y2-3*y3+y4)+
        I*(-20*x1*y1+30*x2*y1-12*x3*y1+x4*y1+30*x1*y2-
        36*x2*y2+9*x3*y2-12*x1*y3+9*x2*y3+x1*y4)+
        F*(-84*x12*y1+168*x1*x2*y1-63*x22*y1-42*x1*x3*y1+
        18*x2*x3*y1+2*x1*x4*y1+84*x12*y2-126*x1*x2*y2+
        27*x22*y2+18*x1*x3*y2-21*x12*y3+18*x1*x2*y3+
        x12*y4)+
        B*(-20*y12+60*y1*y2-36*y22-24*y1*y3+18*y2*y3+
        2*y1*y4)+
        H*(-84*x1*y12+84*x2*y12-21*x3*y12+x4*y12+
        168*x1*y1*y2-126*x2*y1*y2+18*x3*y1*y2-63*x1*y22+
        27*x2*y22-42*x1*y1*y3+18*x2*y1*y3+18*x1*y2*y3+
        2*x1*y1*y4)+
        A*(-84*y13+252*y12*y2-189*y1*y22+27*y23-63*y12*y3+
        54*y1*y2*y3+3*y12*y4);

    var t4 = E*(15*x12-60*x1*x2+54*x22+36*x1*x3-54*x2*x3+
        9*x32-6*x1*x4+6*x2*x4)+
        D*(126*x13-504*x12*x2+567*x1*x22-162*x23+
        189*x12*x3-324*x1*x2*x3+81*x22*x3+27*x1*x32-
        18*x12*x4+18*x1*x2*x4)+
        I*(15*x1*y1-30*x2*y1+18*x3*y1-3*x4*y1-30*x1*y2+
        54*x2*y2-27*x3*y2+3*x4*y2+18*x1*y3-27*x2*y3+
        9*x3*y3-3*x1*y4+3*x2*y4)+
        F*(126*x12*y1-336*x1*x2*y1+189*x22*y1+126*x1*x3*y1-
        108*x2*x3*y1+9*x32*y1-12*x1*x4*y1+6*x2*x4*y1-
        168*x12*y2+378*x1*x2*y2-162*x22*y2-108*x1*x3*y2+
        54*x2*x3*y2+6*x1*x4*y2+63*x12*y3-108*x1*x2*y3+
        27*x22*y3+18*x1*x3*y3-6*x12*y4+6*x1*x2*y4)+
        B*(15*y12-60*y1*y2+54*y22+36*y1*y3-54*y2*y3+9*y32-
        6*y1*y4+6*y2*y4)+
        H*(126*x1*y12-168*x2*y12+63*x3*y12-6*x4*y12-
        336*x1*y1*y2+378*x2*y1*y2-108*x3*y1*y2+6*x4*y1*y2+
        189*x1*y22-162*x2*y22+27*x3*y22+126*x1*y1*y3-
        108*x2*y1*y3+18*x3*y1*y3-108*x1*y2*y3+54*x2*y2*y3+
        9*x1*y32-12*x1*y1*y4+6*x2*y1*y4+6*x1*y2*y4)+
        A*(126*y13-504*y12*y2+567*y1*y22-162*y23+
        189*y12*y3-324*y1*y2*y3+81*y22*y3+27*y1*y32-
        18*y12*y4+18*y1*y2*y4);

    var t5 = E*(-6*x12+30*x1*x2-36*x22-24*x1*x3+54*x2*x3-
        18*x32+6*x1*x4-12*x2*x4+6*x3*x4)+
        D*(-126*x13+630*x12*x2-945*x1*x22+405*x23-
        315*x12*x3+810*x1*x2*x3-405*x22*x3-135*x1*x32+
        81*x2*x32+45*x12*x4-90*x1*x2*x4+27*x22*x4+
        18*x1*x3*x4)+
        I*(-6*x1*y1+15*x2*y1-12*x3*y1+3*x4*y1+15*x1*y2-
        36*x2*y2+27*x3*y2-6*x4*y2-12*x1*y3+27*x2*y3-
        18*x3*y3+3*x4*y3+3*x1*y4-6*x2*y4+3*x3*y4)+
        F*(-126*x12*y1+420*x1*x2*y1-315*x22*y1-210*x1*x3*y1+
        270*x2*x3*y1-45*x32*y1+30*x1*x4*y1-30*x2*x4*y1+
        6*x3*x4*y1+210*x12*y2-630*x1*x2*y2+405*x22*y2+
        270*x1*x3*y2-270*x2*x3*y2+27*x32*y2-30*x1*x4*y2+
        18*x2*x4*y2-105*x12*y3+270*x1*x2*y3-135*x22*y3-
        90*x1*x3*y3+54*x2*x3*y3+6*x1*x4*y3+15*x12*y4-
        30*x1*x2*y4+9*x22*y4+6*x1*x3*y4)+
        B*(-6*y12+30*y1*y2-36*y22-24*y1*y3+54*y2*y3-
        18*y32+6*y1*y4-12*y2*y4+6*y3*y4)+
        H*(-126*x1*y12+210*x2*y12-105*x3*y12+15*x4*y12+
        420*x1*y1*y2-630*x2*y1*y2+270*x3*y1*y2-30*x4*y1*y2-
        315*x1*y22+405*x2*y22-135*x3*y22+9*x4*y22-
        210*x1*y1*y3+270*x2*y1*y3-90*x3*y1*y3+6*x4*y1*y3+
        270*x1*y2*y3-270*x2*y2*y3+54*x3*y2*y3-45*x1*y32+
        27*x2*y32+30*x1*y1*y4-30*x2*y1*y4+6*x3*y1*y4-
        30*x1*y2*y4+18*x2*y2*y4+6*x1*y3*y4)+
        A*(-126*y13+630*y12*y2-945*y1*y22+405*y23-
        315*y12*y3+810*y1*y2*y3-405*y22*y3-135*y1*y32+
        81*y2*y32+45*y12*y4-90*y1*y2*y4+27*y22*y4+
        18*y1*y3*y4);

    var t6 = E*(x12-6*x1*x2+9*x22+6*x1*x3-18*x2*x3+9*x32-
        2*x1*x4+6*x2*x4-6*x3*x4+x42)+
        D*(84*x13-504*x12*x2+945*x1*x22-540*x23+
        315*x12*x3-1080*x1*x2*x3+810*x22*x3+270*x1*x32-
        324*x2*x32+27*x33-60*x12*x4+180*x1*x2*x4-
        108*x22*x4-72*x1*x3*x4+54*x2*x3*x4+3*x1*x42)+
        I*(x1*y1-3*x2*y1+3*x3*y1-x4*y1-3*x1*y2+9*x2*y2-
        9*x3*y2+3*x4*y2+3*x1*y3-9*x2*y3+9*x3*y3-3*x4*y3-
        x1*y4+3*x2*y4-3*x3*y4+x4*y4)+
        F*(84*x12*y1-336*x1*x2*y1+315*x22*y1+210*x1*x3*y1-
        360*x2*x3*y1+90*x32*y1-40*x1*x4*y1+60*x2*x4*y1-
        24*x3*x4*y1+x42*y1-168*x12*y2+630*x1*x2*y2-
        540*x22*y2-360*x1*x3*y2+540*x2*x3*y2-108*x32*y2+
        60*x1*x4*y2-72*x2*x4*y2+18*x3*x4*y2+105*x12*y3-
        360*x1*x2*y3+270*x22*y3+180*x1*x3*y3-216*x2*x3*y3+
        27*x32*y3-24*x1*x4*y3+18*x2*x4*y3-20*x12*y4+
        60*x1*x2*y4-36*x22*y4-24*x1*x3*y4+18*x2*x3*y4+
        2*x1*x4*y4)+
        B*(y12-6*y1*y2+9*y22+6*y1*y3-18*y2*y3+9*y32-
        2*y1*y4+6*y2*y4-6*y3*y4+y42)+
        H*(84*x1*y12-168*x2*y12+105*x3*y12-20*x4*y12-
        336*x1*y1*y2+630*x2*y1*y2-360*x3*y1*y2+60*x4*y1*y2+
        315*x1*y22-540*x2*y22+270*x3*y22-36*x4*y22+
        210*x1*y1*y3-360*x2*y1*y3+180*x3*y1*y3-24*x4*y1*y3-
        360*x1*y2*y3+540*x2*y2*y3-216*x3*y2*y3+18*x4*y2*y3+
        90*x1*y32-108*x2*y32+27*x3*y32-40*x1*y1*y4+
        60*x2*y1*y4-24*x3*y1*y4+2*x4*y1*y4+60*x1*y2*y4-
        72*x2*y2*y4+18*x3*y2*y4-24*x1*y3*y4+18*x2*y3*y4+
        x1*y42)+
        A*(84*y13-504*y12*y2+945*y1*y22-540*y23+
        315*y12*y3-1080*y1*y2*y3+810*y22*y3+270*y1*y32-
        324*y2*y32+27*y33-60*y12*y4+180*y1*y2*y4-
        108*y22*y4-72*y1*y3*y4+54*y2*y3*y4+3*y1*y42);

    var t7 = D*(-36*x13+252*x12*x2-567*x1*x22+405*x23-
        189*x12*x3+810*x1*x2*x3-810*x22*x3-270*x1*x32+
        486*x2*x32-81*x33+45*x12*x4-180*x1*x2*x4+
        162*x22*x4+108*x1*x3*x4-162*x2*x3*x4+27*x32*x4-
        9*x1*x42+9*x2*x42)+
        F*(-36*x12*y1+168*x1*x2*y1-189*x22*y1-126*x1*x3*y1+
        270*x2*x3*y1-90*x32*y1+30*x1*x4*y1-60*x2*x4*y1+
        36*x3*x4*y1-3*x42*y1+84*x12*y2-378*x1*x2*y2+
        405*x22*y2+270*x1*x3*y2-540*x2*x3*y2+162*x32*y2-
        60*x1*x4*y2+108*x2*x4*y2-54*x3*x4*y2+3*x42*y2-
        63*x12*y3+270*x1*x2*y3-270*x22*y3-180*x1*x3*y3+
        324*x2*x3*y3-81*x32*y3+36*x1*x4*y3-54*x2*x4*y3+
        18*x3*x4*y3+15*x12*y4-60*x1*x2*y4+54*x22*y4+
        36*x1*x3*y4-54*x2*x3*y4+9*x32*y4-6*x1*x4*y4+
        6*x2*x4*y4)+
        H*(-36*x1*y12+84*x2*y12-63*x3*y12+15*x4*y12+
        168*x1*y1*y2-378*x2*y1*y2+270*x3*y1*y2-60*x4*y1*y2-
        189*x1*y22+405*x2*y22-270*x3*y22+54*x4*y22-
        126*x1*y1*y3+270*x2*y1*y3-180*x3*y1*y3+36*x4*y1*y3+
        270*x1*y2*y3-540*x2*y2*y3+324*x3*y2*y3-54*x4*y2*y3-
        90*x1*y32+162*x2*y32-81*x3*y32+9*x4*y32+
        30*x1*y1*y4-60*x2*y1*y4+36*x3*y1*y4-6*x4*y1*y4-
        60*x1*y2*y4+108*x2*y2*y4-54*x3*y2*y4+6*x4*y2*y4+
        36*x1*y3*y4-54*x2*y3*y4+18*x3*y3*y4-3*x1*y42+
        3*x2*y42)+
        A*(-36*y13+252*y12*y2-567*y1*y22+405*y23-
        189*y12*y3+810*y1*y2*y3-810*y22*y3-270*y1*y32+
        486*y2*y32-81*y33+45*y12*y4-180*y1*y2*y4+
        162*y22*y4+108*y1*y3*y4-162*y2*y3*y4+27*y32*y4-
        9*y1*y42+9*y2*y42);

    var t8 = D*(9*x13-72*x12*x2+189*x1*x22-162*x23+63*x12*x3-
        324*x1*x2*x3+405*x22*x3+135*x1*x32-324*x2*x32+
        81*x33-18*x12*x4+90*x1*x2*x4-108*x22*x4-
        72*x1*x3*x4+162*x2*x3*x4-54*x32*x4+9*x1*x42-
        18*x2*x42+9*x3*x42)+
        F*(9*x12*y1-48*x1*x2*y1+63*x22*y1+42*x1*x3*y1-
        108*x2*x3*y1+45*x32*y1-12*x1*x4*y1+30*x2*x4*y1-
        24*x3*x4*y1+3*x42*y1-24*x12*y2+126*x1*x2*y2-
        162*x22*y2-108*x1*x3*y2+270*x2*x3*y2-108*x32*y2+
        30*x1*x4*y2-72*x2*x4*y2+54*x3*x4*y2-6*x42*y2+
        21*x12*y3-108*x1*x2*y3+135*x22*y3+90*x1*x3*y3-
        216*x2*x3*y3+81*x32*y3-24*x1*x4*y3+54*x2*x4*y3-
        36*x3*x4*y3+3*x42*y3-6*x12*y4+30*x1*x2*y4-
        36*x22*y4-24*x1*x3*y4+54*x2*x3*y4-18*x32*y4+
        6*x1*x4*y4-12*x2*x4*y4+6*x3*x4*y4)+
        H*(9*x1*y12-24*x2*y12+21*x3*y12-6*x4*y12-
        48*x1*y1*y2+126*x2*y1*y2-108*x3*y1*y2+30*x4*y1*y2+
        63*x1*y22-162*x2*y22+135*x3*y22-36*x4*y22+
        42*x1*y1*y3-108*x2*y1*y3+90*x3*y1*y3-24*x4*y1*y3-
        108*x1*y2*y3+270*x2*y2*y3-216*x3*y2*y3+54*x4*y2*y3+
        45*x1*y32-108*x2*y32+81*x3*y32-18*x4*y32-
        12*x1*y1*y4+30*x2*y1*y4-24*x3*y1*y4+6*x4*y1*y4+
        30*x1*y2*y4-72*x2*y2*y4+54*x3*y2*y4-12*x4*y2*y4-
        24*x1*y3*y4+54*x2*y3*y4-36*x3*y3*y4+6*x4*y3*y4+
        3*x1*y42-6*x2*y42+3*x3*y42)+
        A*(9*y13-72*y12*y2+189*y1*y22-162*y23+63*y12*y3-
        324*y1*y2*y3+405*y22*y3+135*y1*y32-324*y2*y32+
        81*y33-18*y12*y4+90*y1*y2*y4-108*y22*y4-
        72*y1*y3*y4+162*y2*y3*y4-54*y32*y4+9*y1*y42-
        18*y2*y42+9*y3*y42);

    var t9 = D*(-x13+9*x12*x2-27*x1*x22+27*x23-9*x12*x3+
        54*x1*x2*x3-81*x22*x3-27*x1*x32+81*x2*x32-27*x33+
        3*x12*x4-18*x1*x2*x4+27*x22*x4+18*x1*x3*x4-
        54*x2*x3*x4+27*x32*x4-3*x1*x42+9*x2*x42-9*x3*x42+
        x43)+F*(-x12*y1+6*x1*x2*y1-9*x22*y1-6*x1*x3*y1+
        18*x2*x3*y1-9*x32*y1+2*x1*x4*y1-6*x2*x4*y1+
        6*x3*x4*y1-x42*y1+3*x12*y2-18*x1*x2*y2+27*x22*y2+
        18*x1*x3*y2-54*x2*x3*y2+27*x32*y2-6*x1*x4*y2+
        18*x2*x4*y2-18*x3*x4*y2+3*x42*y2-3*x12*y3+
        18*x1*x2*y3-27*x22*y3-18*x1*x3*y3+54*x2*x3*y3-
        27*x32*y3+6*x1*x4*y3-18*x2*x4*y3+18*x3*x4*y3-
        3*x42*y3+x12*y4-6*x1*x2*y4+9*x22*y4+6*x1*x3*y4-
        18*x2*x3*y4+9*x32*y4-2*x1*x4*y4+6*x2*x4*y4-
        6*x3*x4*y4+x42*y4)+
        H*(-x1*y12+3*x2*y12-3*x3*y12+x4*y12+6*x1*y1*y2-
        18*x2*y1*y2+18*x3*y1*y2-6*x4*y1*y2-9*x1*y22+
        27*x2*y22-27*x3*y22+9*x4*y22-6*x1*y1*y3+
        18*x2*y1*y3-18*x3*y1*y3+6*x4*y1*y3+18*x1*y2*y3-
        54*x2*y2*y3+54*x3*y2*y3-18*x4*y2*y3-9*x1*y32+
        27*x2*y32-27*x3*y32+9*x4*y32+2*x1*y1*y4-
        6*x2*y1*y4+6*x3*y1*y4-2*x4*y1*y4-6*x1*y2*y4+
        18*x2*y2*y4-18*x3*y2*y4+6*x4*y2*y4+6*x1*y3*y4-
        18*x2*y3*y4+18*x3*y3*y4-6*x4*y3*y4-x1*y42+3*x2*y42-
        3*x3*y42+x4*y42)+
        A*(-y13+9*y12*y2-27*y1*y22+27*y23-9*y12*y3+
        54*y1*y2*y3-81*y22*y3-27*y1*y32+81*y2*y32-27*y33+
        3*y12*y4-18*y1*y2*y4+27*y22*y4+18*y1*y3*y4-
        54*y2*y3*y4+27*y32*y4-3*y1*y42+9*y2*y42-9*y3*y42+
        y43);

        var min = 1// Math.min( t0, t1, t2, t3, t4, t5, t6, t7, t8, t9 );
        return [ t0/min, t1/min, t2/min, t3/min, t4/min,
                t5/min, t6/min, t7/min, t8/min, t9/min ];
}


/**
 * Horner's scheme using Newton_Raphson to successively
 * find all roots for a n degree polynomial.
 *     http://en.wikipedia.org/wiki/Horner%27s_method#Polynomial_root_finding
 *
 * Newton-Raphson part of this implementation is adapted from
 * paperjs #Numerical.findRoot method
 *
 * @param  {Array} p - Coefficients of polynomial
 * @return {Array}   - Roots in range [0..1]
 */
function findRoots( _p ){
    var n = _p.length-1, i, j;
    var z = 1, fz, dz;
    var roots = [], p = [], dp = [];
    var rootCount = 0, lastSign ,sign, pi, pj;
    lastSign = _p[0] > 0;
    // dp is the derivative of polynomial p w.r.t t
    for (i = 0, j = n; i <= n ; i++, j--){
        pi = _p[i];
        pj = _p[j];
        sign = pj > 0;
        if( i && pj && lastSign !== sign ){
            ++rootCount;
        }
        if( pj ) lastSign = sign;
        p.push( pi );
        dp.push( i * pi );
    }
    // According to Descartes' rule of signs, the number of real positive roots
    // is equal to or less than the change in signs of consecutive non-zero coeff.
    if( !rootCount )
        return roots;
    dp.shift();
    while( n-- ){
        for (i = 0; i <= MAX_ITERATE; i++) {
            fz = evaluateHorner( p, z );
            dz = fz / evaluateHorner( dp, z );
            // Check if we are done
            if (Math.abs(dz) < EPSILON){
                if( z >= 0 && z <= 1 && evaluateHorner(_p, z) < 0.01 ){
                    roots.push( z );
                }
                break;
            }
            // find the next approximation
            z = z - dz;
        }
        // DEBUG: Should we exit if the root fails to converge?
        if( i > MAX_ITERATE ) break;
        // Divide p itself instead of creating a new array every time. Much faster this way!
        polyLongDivideSelf( p, z );
        dp.pop();
        for (i = 1; i <= n; i++)
            dp[i-1] = i * p[i];
    }
    return roots;
}


/**
*   Adapted From http://www.kevlindev.com/gui/math/intersection/index.htm
*   copyright 2002-2003, Kevin Lindsey
*/
var intersectCurveCurve = function(v1, v2, curve1, curve2, locations) {
    var ax, bx, cx, dx, ay, by, cy, dy;
    var c13x, c12x, c11x, c10x, c13y, c12y, c11y, c10y;
    var c23x, c22x, c21x, c20x, c23y, c22y, c21y, c20y;
    var a1x = v1[0], a1y = v1[1], a2x = v1[2], a2y = v1[3];
    var a3x = v1[4], a3y = v1[5], a4x = v1[6], a4y = v1[7];
    var b1x = v2[0], b1y = v2[1], b2x = v2[2], b2y = v2[3];
    var b3x = v2[4], b3y = v2[5], b4x = v2[6], b4y = v2[7];
    // Calculate the coefficients of cubic polynomial
    ax = a1x*-1; ay = a1y*-1;
    bx = a2x*3; by = a2y*3;
    cx = a3x*-3; cy = a3y*-3;
    c13x = ax+bx+cx+a4x; c13y = ay+by+cy+a4y;
    ax = a1x*3; ay = a1y*3;
    bx = a2x*-6; by = a2y*-6;
    cx = a3x*3; cy = a3y*3;
    c12x = ax+bx+cx; c12y = ay+by+cy;
    ax = a1x*-3; ay = a1y*-3;
    bx = a2x*3; by = a2y*3;
    c11x = ax+bx; c11y = ay+by;
    c10x = a1x; c10y = a1y;
    ax = b1x*-1; ay = b1y*-1;
    bx = b2x*3; by = b2y*3;
    cx = b3x*-3; cy = b3y*-3;
    c23x = ax+bx+cx+b4x; c23y = ay+by+cy+b4y;
    ax = b1x*3; ay = b1y*3;
    bx = b2x*-6; by = b2y*-6;
    cx = b3x*3; cy = b3y*3;
    c22x = ax+bx+cx; c22y = ay+by+cy;
    ax = b1x*-3; ay = b1y*-3;
    bx = b2x*3; by = b2y*3;
    c21x =  ax + bx; c21y =  ay + by;
    c20x = b1x; c20y = b1y;

    var c10x2 = c10x*c10x;
    var c10x3 = c10x*c10x*c10x;
    var c10y2 = c10y*c10y;
    var c10y3 = c10y*c10y*c10y;
    var c11x2 = c11x*c11x;
    var c11x3 = c11x*c11x*c11x;
    var c11y2 = c11y*c11y;
    var c11y3 = c11y*c11y*c11y;
    var c12x2 = c12x*c12x;
    var c12x3 = c12x*c12x*c12x;
    var c12y2 = c12y*c12y;
    var c12y3 = c12y*c12y*c12y;
    var c13x2 = c13x*c13x;
    var c13x3 = c13x*c13x*c13x;
    var c13y2 = c13y*c13y;
    var c13y3 = c13y*c13y*c13y;
    var c20x2 = c20x*c20x;
    var c20x3 = c20x*c20x*c20x;
    var c20y2 = c20y*c20y;
    var c20y3 = c20y*c20y*c20y;
    var c21x2 = c21x*c21x;
    var c21x3 = c21x*c21x*c21x;
    var c21y2 = c21y*c21y;
    var c22x2 = c22x*c22x;
    var c22x3 = c22x*c22x*c22x;
    var c22y2 = c22y*c22y;
    var c23x2 = c23x*c23x;
    var c23x3 = c23x*c23x*c23x;
    var c23y2 = c23y*c23y;
    var c23y3 = c23y*c23y*c23y;
    var poly = new Polynomial(
        -c13x3*c23y3 + c13y3*c23x3 - 3*c13x*c13y2*c23x2*c23y +
            3*c13x2*c13y*c23x*c23y2,
        -6*c13x*c22x*c13y2*c23x*c23y + 6*c13x2*c13y*c22y*c23x*c23y + 3*c22x*c13y3*c23x2 -
            3*c13x3*c22y*c23y2 - 3*c13x*c13y2*c22y*c23x2 + 3*c13x2*c22x*c13y*c23y2,
        -6*c21x*c13x*c13y2*c23x*c23y - 6*c13x*c22x*c13y2*c22y*c23x + 6*c13x2*c22x*c13y*c22y*c23y +
            3*c21x*c13y3*c23x2 + 3*c22x2*c13y3*c23x + 3*c21x*c13x2*c13y*c23y2 - 3*c13x*c21y*c13y2*c23x2 -
            3*c13x*c22x2*c13y2*c23y + c13x2*c13y*c23x*(6*c21y*c23y + 3*c22y2) + c13x3*(-c21y*c23y2 -
            2*c22y2*c23y - c23y*(2*c21y*c23y + c22y2)),
        c11x*c12y*c13x*c13y*c23x*c23y - c11y*c12x*c13x*c13y*c23x*c23y + 6*c21x*c22x*c13y3*c23x +
            3*c11x*c12x*c13x*c13y*c23y2 + 6*c10x*c13x*c13y2*c23x*c23y - 3*c11x*c12x*c13y2*c23x*c23y -
            3*c11y*c12y*c13x*c13y*c23x2 - 6*c10y*c13x2*c13y*c23x*c23y - 6*c20x*c13x*c13y2*c23x*c23y +
            3*c11y*c12y*c13x2*c23x*c23y - 2*c12x*c12y2*c13x*c23x*c23y - 6*c21x*c13x*c22x*c13y2*c23y -
            6*c21x*c13x*c13y2*c22y*c23x - 6*c13x*c21y*c22x*c13y2*c23x + 6*c21x*c13x2*c13y*c22y*c23y +
            2*c12x2*c12y*c13y*c23x*c23y + c22x3*c13y3 - 3*c10x*c13y3*c23x2 + 3*c10y*c13x3*c23y2 +
            3*c20x*c13y3*c23x2 + c12y3*c13x*c23x2 - c12x3*c13y*c23y2 - 3*c10x*c13x2*c13y*c23y2 +
            3*c10y*c13x*c13y2*c23x2 - 2*c11x*c12y*c13x2*c23y2 + c11x*c12y*c13y2*c23x2 - c11y*c12x*c13x2*c23y2 +
            2*c11y*c12x*c13y2*c23x2 + 3*c20x*c13x2*c13y*c23y2 - c12x*c12y2*c13y*c23x2 -
            3*c20y*c13x*c13y2*c23x2 + c12x2*c12y*c13x*c23y2 - 3*c13x*c22x2*c13y2*c22y +
            c13x2*c13y*c23x*(6*c20y*c23y + 6*c21y*c22y) + c13x2*c22x*c13y*(6*c21y*c23y + 3*c22y2) +
            c13x3*(-2*c21y*c22y*c23y - c20y*c23y2 - c22y*(2*c21y*c23y + c22y2) - c23y*(2*c20y*c23y + 2*c21y*c22y)),
        6*c11x*c12x*c13x*c13y*c22y*c23y + c11x*c12y*c13x*c22x*c13y*c23y + c11x*c12y*c13x*c13y*c22y*c23x -
            c11y*c12x*c13x*c22x*c13y*c23y - c11y*c12x*c13x*c13y*c22y*c23x - 6*c11y*c12y*c13x*c22x*c13y*c23x -
            6*c10x*c22x*c13y3*c23x + 6*c20x*c22x*c13y3*c23x + 6*c10y*c13x3*c22y*c23y + 2*c12y3*c13x*c22x*c23x -
            2*c12x3*c13y*c22y*c23y + 6*c10x*c13x*c22x*c13y2*c23y + 6*c10x*c13x*c13y2*c22y*c23x +
            6*c10y*c13x*c22x*c13y2*c23x - 3*c11x*c12x*c22x*c13y2*c23y - 3*c11x*c12x*c13y2*c22y*c23x +
            2*c11x*c12y*c22x*c13y2*c23x + 4*c11y*c12x*c22x*c13y2*c23x - 6*c10x*c13x2*c13y*c22y*c23y -
            6*c10y*c13x2*c22x*c13y*c23y - 6*c10y*c13x2*c13y*c22y*c23x - 4*c11x*c12y*c13x2*c22y*c23y -
            6*c20x*c13x*c22x*c13y2*c23y - 6*c20x*c13x*c13y2*c22y*c23x - 2*c11y*c12x*c13x2*c22y*c23y +
            3*c11y*c12y*c13x2*c22x*c23y + 3*c11y*c12y*c13x2*c22y*c23x - 2*c12x*c12y2*c13x*c22x*c23y -
            2*c12x*c12y2*c13x*c22y*c23x - 2*c12x*c12y2*c22x*c13y*c23x - 6*c20y*c13x*c22x*c13y2*c23x -
            6*c21x*c13x*c21y*c13y2*c23x - 6*c21x*c13x*c22x*c13y2*c22y + 6*c20x*c13x2*c13y*c22y*c23y +
            2*c12x2*c12y*c13x*c22y*c23y + 2*c12x2*c12y*c22x*c13y*c23y + 2*c12x2*c12y*c13y*c22y*c23x +
            3*c21x*c22x2*c13y3 + 3*c21x2*c13y3*c23x - 3*c13x*c21y*c22x2*c13y2 - 3*c21x2*c13x*c13y2*c23y +
            c13x2*c22x*c13y*(6*c20y*c23y + 6*c21y*c22y) + c13x2*c13y*c23x*(6*c20y*c22y + 3*c21y2) +
            c21x*c13x2*c13y*(6*c21y*c23y + 3*c22y2) + c13x3*(-2*c20y*c22y*c23y - c23y*(2*c20y*c22y + c21y2) -
            c21y*(2*c21y*c23y + c22y2) - c22y*(2*c20y*c23y + 2*c21y*c22y)),
        c11x*c21x*c12y*c13x*c13y*c23y + c11x*c12y*c13x*c21y*c13y*c23x + c11x*c12y*c13x*c22x*c13y*c22y -
            c11y*c12x*c21x*c13x*c13y*c23y - c11y*c12x*c13x*c21y*c13y*c23x - c11y*c12x*c13x*c22x*c13y*c22y -
            6*c11y*c21x*c12y*c13x*c13y*c23x - 6*c10x*c21x*c13y3*c23x + 6*c20x*c21x*c13y3*c23x +
            2*c21x*c12y3*c13x*c23x + 6*c10x*c21x*c13x*c13y2*c23y + 6*c10x*c13x*c21y*c13y2*c23x +
            6*c10x*c13x*c22x*c13y2*c22y + 6*c10y*c21x*c13x*c13y2*c23x - 3*c11x*c12x*c21x*c13y2*c23y -
            3*c11x*c12x*c21y*c13y2*c23x - 3*c11x*c12x*c22x*c13y2*c22y + 2*c11x*c21x*c12y*c13y2*c23x +
            4*c11y*c12x*c21x*c13y2*c23x - 6*c10y*c21x*c13x2*c13y*c23y - 6*c10y*c13x2*c21y*c13y*c23x -
            6*c10y*c13x2*c22x*c13y*c22y - 6*c20x*c21x*c13x*c13y2*c23y - 6*c20x*c13x*c21y*c13y2*c23x -
            6*c20x*c13x*c22x*c13y2*c22y + 3*c11y*c21x*c12y*c13x2*c23y - 3*c11y*c12y*c13x*c22x2*c13y +
            3*c11y*c12y*c13x2*c21y*c23x + 3*c11y*c12y*c13x2*c22x*c22y - 2*c12x*c21x*c12y2*c13x*c23y -
            2*c12x*c21x*c12y2*c13y*c23x - 2*c12x*c12y2*c13x*c21y*c23x - 2*c12x*c12y2*c13x*c22x*c22y -
            6*c20y*c21x*c13x*c13y2*c23x - 6*c21x*c13x*c21y*c22x*c13y2 + 6*c20y*c13x2*c21y*c13y*c23x +
            2*c12x2*c21x*c12y*c13y*c23y + 2*c12x2*c12y*c21y*c13y*c23x + 2*c12x2*c12y*c22x*c13y*c22y -
            3*c10x*c22x2*c13y3 + 3*c20x*c22x2*c13y3 + 3*c21x2*c22x*c13y3 + c12y3*c13x*c22x2 +
            3*c10y*c13x*c22x2*c13y2 + c11x*c12y*c22x2*c13y2 + 2*c11y*c12x*c22x2*c13y2 -
            c12x*c12y2*c22x2*c13y - 3*c20y*c13x*c22x2*c13y2 - 3*c21x2*c13x*c13y2*c22y +
            c12x2*c12y*c13x*(2*c21y*c23y + c22y2) + c11x*c12x*c13x*c13y*(6*c21y*c23y + 3*c22y2) +
            c21x*c13x2*c13y*(6*c20y*c23y + 6*c21y*c22y) + c12x3*c13y*(-2*c21y*c23y - c22y2) +
            c10y*c13x3*(6*c21y*c23y + 3*c22y2) + c11y*c12x*c13x2*(-2*c21y*c23y - c22y2) +
            c11x*c12y*c13x2*(-4*c21y*c23y - 2*c22y2) + c10x*c13x2*c13y*(-6*c21y*c23y - 3*c22y2) +
            c13x2*c22x*c13y*(6*c20y*c22y + 3*c21y2) + c20x*c13x2*c13y*(6*c21y*c23y + 3*c22y2) +
            c13x3*(-2*c20y*c21y*c23y - c22y*(2*c20y*c22y + c21y2) - c20y*(2*c21y*c23y + c22y2) -
            c21y*(2*c20y*c23y + 2*c21y*c22y)),
        -c10x*c11x*c12y*c13x*c13y*c23y + c10x*c11y*c12x*c13x*c13y*c23y + 6*c10x*c11y*c12y*c13x*c13y*c23x -
            6*c10y*c11x*c12x*c13x*c13y*c23y - c10y*c11x*c12y*c13x*c13y*c23x + c10y*c11y*c12x*c13x*c13y*c23x +
            c11x*c11y*c12x*c12y*c13x*c23y - c11x*c11y*c12x*c12y*c13y*c23x + c11x*c20x*c12y*c13x*c13y*c23y +
            c11x*c20y*c12y*c13x*c13y*c23x + c11x*c21x*c12y*c13x*c13y*c22y + c11x*c12y*c13x*c21y*c22x*c13y -
            c20x*c11y*c12x*c13x*c13y*c23y - 6*c20x*c11y*c12y*c13x*c13y*c23x - c11y*c12x*c20y*c13x*c13y*c23x -
            c11y*c12x*c21x*c13x*c13y*c22y - c11y*c12x*c13x*c21y*c22x*c13y - 6*c11y*c21x*c12y*c13x*c22x*c13y -
            6*c10x*c20x*c13y3*c23x - 6*c10x*c21x*c22x*c13y3 - 2*c10x*c12y3*c13x*c23x + 6*c20x*c21x*c22x*c13y3 +
            2*c20x*c12y3*c13x*c23x + 2*c21x*c12y3*c13x*c22x + 2*c10y*c12x3*c13y*c23y - 6*c10x*c10y*c13x*c13y2*c23x +
            3*c10x*c11x*c12x*c13y2*c23y - 2*c10x*c11x*c12y*c13y2*c23x - 4*c10x*c11y*c12x*c13y2*c23x +
            3*c10y*c11x*c12x*c13y2*c23x + 6*c10x*c10y*c13x2*c13y*c23y + 6*c10x*c20x*c13x*c13y2*c23y -
            3*c10x*c11y*c12y*c13x2*c23y + 2*c10x*c12x*c12y2*c13x*c23y + 2*c10x*c12x*c12y2*c13y*c23x +
            6*c10x*c20y*c13x*c13y2*c23x + 6*c10x*c21x*c13x*c13y2*c22y + 6*c10x*c13x*c21y*c22x*c13y2 +
            4*c10y*c11x*c12y*c13x2*c23y + 6*c10y*c20x*c13x*c13y2*c23x + 2*c10y*c11y*c12x*c13x2*c23y -
            3*c10y*c11y*c12y*c13x2*c23x + 2*c10y*c12x*c12y2*c13x*c23x + 6*c10y*c21x*c13x*c22x*c13y2 -
            3*c11x*c20x*c12x*c13y2*c23y + 2*c11x*c20x*c12y*c13y2*c23x + c11x*c11y*c12y2*c13x*c23x -
            3*c11x*c12x*c20y*c13y2*c23x - 3*c11x*c12x*c21x*c13y2*c22y - 3*c11x*c12x*c21y*c22x*c13y2 +
            2*c11x*c21x*c12y*c22x*c13y2 + 4*c20x*c11y*c12x*c13y2*c23x + 4*c11y*c12x*c21x*c22x*c13y2 -
            2*c10x*c12x2*c12y*c13y*c23y - 6*c10y*c20x*c13x2*c13y*c23y - 6*c10y*c20y*c13x2*c13y*c23x -
            6*c10y*c21x*c13x2*c13y*c22y - 2*c10y*c12x2*c12y*c13x*c23y - 2*c10y*c12x2*c12y*c13y*c23x -
            6*c10y*c13x2*c21y*c22x*c13y - c11x*c11y*c12x2*c13y*c23y - 2*c11x*c11y2*c13x*c13y*c23x +
            3*c20x*c11y*c12y*c13x2*c23y - 2*c20x*c12x*c12y2*c13x*c23y - 2*c20x*c12x*c12y2*c13y*c23x -
            6*c20x*c20y*c13x*c13y2*c23x - 6*c20x*c21x*c13x*c13y2*c22y - 6*c20x*c13x*c21y*c22x*c13y2 +
            3*c11y*c20y*c12y*c13x2*c23x + 3*c11y*c21x*c12y*c13x2*c22y + 3*c11y*c12y*c13x2*c21y*c22x -
            2*c12x*c20y*c12y2*c13x*c23x - 2*c12x*c21x*c12y2*c13x*c22y - 2*c12x*c21x*c12y2*c22x*c13y -
            2*c12x*c12y2*c13x*c21y*c22x - 6*c20y*c21x*c13x*c22x*c13y2 - c11y2*c12x*c12y*c13x*c23x +
            2*c20x*c12x2*c12y*c13y*c23y + 6*c20y*c13x2*c21y*c22x*c13y + 2*c11x2*c11y*c13x*c13y*c23y +
            c11x2*c12x*c12y*c13y*c23y + 2*c12x2*c20y*c12y*c13y*c23x + 2*c12x2*c21x*c12y*c13y*c22y +
            2*c12x2*c12y*c21y*c22x*c13y + c21x3*c13y3 + 3*c10x2*c13y3*c23x - 3*c10y2*c13x3*c23y +
            3*c20x2*c13y3*c23x + c11y3*c13x2*c23x - c11x3*c13y2*c23y - c11x*c11y2*c13x2*c23y +
            c11x2*c11y*c13y2*c23x - 3*c10x2*c13x*c13y2*c23y + 3*c10y2*c13x2*c13y*c23x - c11x2*c12y2*c13x*c23y +
            c11y2*c12x2*c13y*c23x - 3*c21x2*c13x*c21y*c13y2 - 3*c20x2*c13x*c13y2*c23y + 3*c20y2*c13x2*c13y*c23x +
            c11x*c12x*c13x*c13y*(6*c20y*c23y + 6*c21y*c22y) + c12x3*c13y*(-2*c20y*c23y - 2*c21y*c22y) +
            c10y*c13x3*(6*c20y*c23y + 6*c21y*c22y) + c11y*c12x*c13x2*(-2*c20y*c23y - 2*c21y*c22y) +
            c12x2*c12y*c13x*(2*c20y*c23y + 2*c21y*c22y) + c11x*c12y*c13x2*(-4*c20y*c23y - 4*c21y*c22y) +
            c10x*c13x2*c13y*(-6*c20y*c23y - 6*c21y*c22y) + c20x*c13x2*c13y*(6*c20y*c23y + 6*c21y*c22y) +
            c21x*c13x2*c13y*(6*c20y*c22y + 3*c21y2) + c13x3*(-2*c20y*c21y*c22y - c20y2*c23y -
            c21y*(2*c20y*c22y + c21y2) - c20y*(2*c20y*c23y + 2*c21y*c22y)),
        -c10x*c11x*c12y*c13x*c13y*c22y + c10x*c11y*c12x*c13x*c13y*c22y + 6*c10x*c11y*c12y*c13x*c22x*c13y -
            6*c10y*c11x*c12x*c13x*c13y*c22y - c10y*c11x*c12y*c13x*c22x*c13y + c10y*c11y*c12x*c13x*c22x*c13y +
            c11x*c11y*c12x*c12y*c13x*c22y - c11x*c11y*c12x*c12y*c22x*c13y + c11x*c20x*c12y*c13x*c13y*c22y +
            c11x*c20y*c12y*c13x*c22x*c13y + c11x*c21x*c12y*c13x*c21y*c13y - c20x*c11y*c12x*c13x*c13y*c22y -
            6*c20x*c11y*c12y*c13x*c22x*c13y - c11y*c12x*c20y*c13x*c22x*c13y - c11y*c12x*c21x*c13x*c21y*c13y -
            6*c10x*c20x*c22x*c13y3 - 2*c10x*c12y3*c13x*c22x + 2*c20x*c12y3*c13x*c22x + 2*c10y*c12x3*c13y*c22y -
            6*c10x*c10y*c13x*c22x*c13y2 + 3*c10x*c11x*c12x*c13y2*c22y - 2*c10x*c11x*c12y*c22x*c13y2 -
            4*c10x*c11y*c12x*c22x*c13y2 + 3*c10y*c11x*c12x*c22x*c13y2 + 6*c10x*c10y*c13x2*c13y*c22y +
            6*c10x*c20x*c13x*c13y2*c22y - 3*c10x*c11y*c12y*c13x2*c22y + 2*c10x*c12x*c12y2*c13x*c22y +
            2*c10x*c12x*c12y2*c22x*c13y + 6*c10x*c20y*c13x*c22x*c13y2 + 6*c10x*c21x*c13x*c21y*c13y2 +
            4*c10y*c11x*c12y*c13x2*c22y + 6*c10y*c20x*c13x*c22x*c13y2 + 2*c10y*c11y*c12x*c13x2*c22y -
            3*c10y*c11y*c12y*c13x2*c22x + 2*c10y*c12x*c12y2*c13x*c22x - 3*c11x*c20x*c12x*c13y2*c22y +
            2*c11x*c20x*c12y*c22x*c13y2 + c11x*c11y*c12y2*c13x*c22x - 3*c11x*c12x*c20y*c22x*c13y2 -
            3*c11x*c12x*c21x*c21y*c13y2 + 4*c20x*c11y*c12x*c22x*c13y2 - 2*c10x*c12x2*c12y*c13y*c22y -
            6*c10y*c20x*c13x2*c13y*c22y - 6*c10y*c20y*c13x2*c22x*c13y - 6*c10y*c21x*c13x2*c21y*c13y -
            2*c10y*c12x2*c12y*c13x*c22y - 2*c10y*c12x2*c12y*c22x*c13y - c11x*c11y*c12x2*c13y*c22y -
            2*c11x*c11y2*c13x*c22x*c13y + 3*c20x*c11y*c12y*c13x2*c22y - 2*c20x*c12x*c12y2*c13x*c22y -
            2*c20x*c12x*c12y2*c22x*c13y - 6*c20x*c20y*c13x*c22x*c13y2 - 6*c20x*c21x*c13x*c21y*c13y2 +
            3*c11y*c20y*c12y*c13x2*c22x + 3*c11y*c21x*c12y*c13x2*c21y - 2*c12x*c20y*c12y2*c13x*c22x -
            2*c12x*c21x*c12y2*c13x*c21y - c11y2*c12x*c12y*c13x*c22x + 2*c20x*c12x2*c12y*c13y*c22y -
            3*c11y*c21x2*c12y*c13x*c13y + 6*c20y*c21x*c13x2*c21y*c13y + 2*c11x2*c11y*c13x*c13y*c22y +
            c11x2*c12x*c12y*c13y*c22y + 2*c12x2*c20y*c12y*c22x*c13y + 2*c12x2*c21x*c12y*c21y*c13y -
            3*c10x*c21x2*c13y3 + 3*c20x*c21x2*c13y3 + 3*c10x2*c22x*c13y3 - 3*c10y2*c13x3*c22y + 3*c20x2*c22x*c13y3 +
            c21x2*c12y3*c13x + c11y3*c13x2*c22x - c11x3*c13y2*c22y + 3*c10y*c21x2*c13x*c13y2 -
            c11x*c11y2*c13x2*c22y + c11x*c21x2*c12y*c13y2 + 2*c11y*c12x*c21x2*c13y2 + c11x2*c11y*c22x*c13y2 -
            c12x*c21x2*c12y2*c13y - 3*c20y*c21x2*c13x*c13y2 - 3*c10x2*c13x*c13y2*c22y + 3*c10y2*c13x2*c22x*c13y -
            c11x2*c12y2*c13x*c22y + c11y2*c12x2*c22x*c13y - 3*c20x2*c13x*c13y2*c22y + 3*c20y2*c13x2*c22x*c13y +
            c12x2*c12y*c13x*(2*c20y*c22y + c21y2) + c11x*c12x*c13x*c13y*(6*c20y*c22y + 3*c21y2) +
            c12x3*c13y*(-2*c20y*c22y - c21y2) + c10y*c13x3*(6*c20y*c22y + 3*c21y2) +
            c11y*c12x*c13x2*(-2*c20y*c22y - c21y2) + c11x*c12y*c13x2*(-4*c20y*c22y - 2*c21y2) +
            c10x*c13x2*c13y*(-6*c20y*c22y - 3*c21y2) + c20x*c13x2*c13y*(6*c20y*c22y + 3*c21y2) +
            c13x3*(-2*c20y*c21y2 - c20y2*c22y - c20y*(2*c20y*c22y + c21y2)),
        -c10x*c11x*c12y*c13x*c21y*c13y + c10x*c11y*c12x*c13x*c21y*c13y + 6*c10x*c11y*c21x*c12y*c13x*c13y -
            6*c10y*c11x*c12x*c13x*c21y*c13y - c10y*c11x*c21x*c12y*c13x*c13y + c10y*c11y*c12x*c21x*c13x*c13y -
            c11x*c11y*c12x*c21x*c12y*c13y + c11x*c11y*c12x*c12y*c13x*c21y + c11x*c20x*c12y*c13x*c21y*c13y +
            6*c11x*c12x*c20y*c13x*c21y*c13y + c11x*c20y*c21x*c12y*c13x*c13y - c20x*c11y*c12x*c13x*c21y*c13y -
            6*c20x*c11y*c21x*c12y*c13x*c13y - c11y*c12x*c20y*c21x*c13x*c13y - 6*c10x*c20x*c21x*c13y3 -
            2*c10x*c21x*c12y3*c13x + 6*c10y*c20y*c13x3*c21y + 2*c20x*c21x*c12y3*c13x + 2*c10y*c12x3*c21y*c13y -
            2*c12x3*c20y*c21y*c13y - 6*c10x*c10y*c21x*c13x*c13y2 + 3*c10x*c11x*c12x*c21y*c13y2 -
            2*c10x*c11x*c21x*c12y*c13y2 - 4*c10x*c11y*c12x*c21x*c13y2 + 3*c10y*c11x*c12x*c21x*c13y2 +
            6*c10x*c10y*c13x2*c21y*c13y + 6*c10x*c20x*c13x*c21y*c13y2 - 3*c10x*c11y*c12y*c13x2*c21y +
            2*c10x*c12x*c21x*c12y2*c13y + 2*c10x*c12x*c12y2*c13x*c21y + 6*c10x*c20y*c21x*c13x*c13y2 +
            4*c10y*c11x*c12y*c13x2*c21y + 6*c10y*c20x*c21x*c13x*c13y2 + 2*c10y*c11y*c12x*c13x2*c21y -
            3*c10y*c11y*c21x*c12y*c13x2 + 2*c10y*c12x*c21x*c12y2*c13x - 3*c11x*c20x*c12x*c21y*c13y2 +
            2*c11x*c20x*c21x*c12y*c13y2 + c11x*c11y*c21x*c12y2*c13x - 3*c11x*c12x*c20y*c21x*c13y2 +
            4*c20x*c11y*c12x*c21x*c13y2 - 6*c10x*c20y*c13x2*c21y*c13y - 2*c10x*c12x2*c12y*c21y*c13y -
            6*c10y*c20x*c13x2*c21y*c13y - 6*c10y*c20y*c21x*c13x2*c13y - 2*c10y*c12x2*c21x*c12y*c13y -
            2*c10y*c12x2*c12y*c13x*c21y - c11x*c11y*c12x2*c21y*c13y - 4*c11x*c20y*c12y*c13x2*c21y -
            2*c11x*c11y2*c21x*c13x*c13y + 3*c20x*c11y*c12y*c13x2*c21y - 2*c20x*c12x*c21x*c12y2*c13y -
            2*c20x*c12x*c12y2*c13x*c21y - 6*c20x*c20y*c21x*c13x*c13y2 - 2*c11y*c12x*c20y*c13x2*c21y +
            3*c11y*c20y*c21x*c12y*c13x2 - 2*c12x*c20y*c21x*c12y2*c13x - c11y2*c12x*c21x*c12y*c13x +
            6*c20x*c20y*c13x2*c21y*c13y + 2*c20x*c12x2*c12y*c21y*c13y + 2*c11x2*c11y*c13x*c21y*c13y +
            c11x2*c12x*c12y*c21y*c13y + 2*c12x2*c20y*c21x*c12y*c13y + 2*c12x2*c20y*c12y*c13x*c21y +
            3*c10x2*c21x*c13y3 - 3*c10y2*c13x3*c21y + 3*c20x2*c21x*c13y3 + c11y3*c21x*c13x2 - c11x3*c21y*c13y2 -
            3*c20y2*c13x3*c21y - c11x*c11y2*c13x2*c21y + c11x2*c11y*c21x*c13y2 - 3*c10x2*c13x*c21y*c13y2 +
            3*c10y2*c21x*c13x2*c13y - c11x2*c12y2*c13x*c21y + c11y2*c12x2*c21x*c13y - 3*c20x2*c13x*c21y*c13y2 +
            3*c20y2*c21x*c13x2*c13y,
        c10x*c10y*c11x*c12y*c13x*c13y - c10x*c10y*c11y*c12x*c13x*c13y + c10x*c11x*c11y*c12x*c12y*c13y -
            c10y*c11x*c11y*c12x*c12y*c13x - c10x*c11x*c20y*c12y*c13x*c13y + 6*c10x*c20x*c11y*c12y*c13x*c13y +
            c10x*c11y*c12x*c20y*c13x*c13y - c10y*c11x*c20x*c12y*c13x*c13y - 6*c10y*c11x*c12x*c20y*c13x*c13y +
            c10y*c20x*c11y*c12x*c13x*c13y - c11x*c20x*c11y*c12x*c12y*c13y + c11x*c11y*c12x*c20y*c12y*c13x +
            c11x*c20x*c20y*c12y*c13x*c13y - c20x*c11y*c12x*c20y*c13x*c13y - 2*c10x*c20x*c12y3*c13x +
            2*c10y*c12x3*c20y*c13y - 3*c10x*c10y*c11x*c12x*c13y2 - 6*c10x*c10y*c20x*c13x*c13y2 +
            3*c10x*c10y*c11y*c12y*c13x2 - 2*c10x*c10y*c12x*c12y2*c13x - 2*c10x*c11x*c20x*c12y*c13y2 -
            c10x*c11x*c11y*c12y2*c13x + 3*c10x*c11x*c12x*c20y*c13y2 - 4*c10x*c20x*c11y*c12x*c13y2 +
            3*c10y*c11x*c20x*c12x*c13y2 + 6*c10x*c10y*c20y*c13x2*c13y + 2*c10x*c10y*c12x2*c12y*c13y +
            2*c10x*c11x*c11y2*c13x*c13y + 2*c10x*c20x*c12x*c12y2*c13y + 6*c10x*c20x*c20y*c13x*c13y2 -
            3*c10x*c11y*c20y*c12y*c13x2 + 2*c10x*c12x*c20y*c12y2*c13x + c10x*c11y2*c12x*c12y*c13x +
            c10y*c11x*c11y*c12x2*c13y + 4*c10y*c11x*c20y*c12y*c13x2 - 3*c10y*c20x*c11y*c12y*c13x2 +
            2*c10y*c20x*c12x*c12y2*c13x + 2*c10y*c11y*c12x*c20y*c13x2 + c11x*c20x*c11y*c12y2*c13x -
            3*c11x*c20x*c12x*c20y*c13y2 - 2*c10x*c12x2*c20y*c12y*c13y - 6*c10y*c20x*c20y*c13x2*c13y -
            2*c10y*c20x*c12x2*c12y*c13y - 2*c10y*c11x2*c11y*c13x*c13y - c10y*c11x2*c12x*c12y*c13y -
            2*c10y*c12x2*c20y*c12y*c13x - 2*c11x*c20x*c11y2*c13x*c13y - c11x*c11y*c12x2*c20y*c13y +
            3*c20x*c11y*c20y*c12y*c13x2 - 2*c20x*c12x*c20y*c12y2*c13x - c20x*c11y2*c12x*c12y*c13x +
            3*c10y2*c11x*c12x*c13x*c13y + 3*c11x*c12x*c20y2*c13x*c13y + 2*c20x*c12x2*c20y*c12y*c13y -
            3*c10x2*c11y*c12y*c13x*c13y + 2*c11x2*c11y*c20y*c13x*c13y + c11x2*c12x*c20y*c12y*c13y -
            3*c20x2*c11y*c12y*c13x*c13y - c10x3*c13y3 + c10y3*c13x3 + c20x3*c13y3 - c20y3*c13x3 -
            3*c10x*c20x2*c13y3 - c10x*c11y3*c13x2 + 3*c10x2*c20x*c13y3 + c10y*c11x3*c13y2 +
            3*c10y*c20y2*c13x3 + c20x*c11y3*c13x2 + c10x2*c12y3*c13x - 3*c10y2*c20y*c13x3 - c10y2*c12x3*c13y +
            c20x2*c12y3*c13x - c11x3*c20y*c13y2 - c12x3*c20y2*c13y - c10x*c11x2*c11y*c13y2 +
            c10y*c11x*c11y2*c13x2 - 3*c10x*c10y2*c13x2*c13y - c10x*c11y2*c12x2*c13y + c10y*c11x2*c12y2*c13x -
            c11x*c11y2*c20y*c13x2 + 3*c10x2*c10y*c13x*c13y2 + c10x2*c11x*c12y*c13y2 +
            2*c10x2*c11y*c12x*c13y2 - 2*c10y2*c11x*c12y*c13x2 - c10y2*c11y*c12x*c13x2 + c11x2*c20x*c11y*c13y2 -
            3*c10x*c20y2*c13x2*c13y + 3*c10y*c20x2*c13x*c13y2 + c11x*c20x2*c12y*c13y2 - 2*c11x*c20y2*c12y*c13x2 +
            c20x*c11y2*c12x2*c13y - c11y*c12x*c20y2*c13x2 - c10x2*c12x*c12y2*c13y - 3*c10x2*c20y*c13x*c13y2 +
            3*c10y2*c20x*c13x2*c13y + c10y2*c12x2*c12y*c13x - c11x2*c20y*c12y2*c13x + 2*c20x2*c11y*c12x*c13y2 +
            3*c20x*c20y2*c13x2*c13y - c20x2*c12x*c12y2*c13y - 3*c20x2*c20y*c13x*c13y2 + c12x2*c20y2*c12y*c13x
    );

    var roots = poly.getRootsInInterval(0,1);
    for ( var i = 0; i < roots.length; i++ ) {
        var s = roots[i];
        var xRoots = new Polynomial(
            c13x,
            c12x,
            c11x,
            c10x - c20x - s*c21x - s*s*c22x - s*s*s*c23x
        ).getRoots();
        var yRoots = new Polynomial(
            c13y,
            c12y,
            c11y,
            c10y - c20y - s*c21y - s*s*c22y - s*s*s*c23y
        ).getRoots();

        if ( xRoots.length > 0 && yRoots.length > 0 ) {
            var TOLERANCE = 1e-4;

            checkRoots:
            for ( var j = 0; j < xRoots.length; j++ ) {
                var xRoot = xRoots[j];

                if ( 0 <= xRoot && xRoot <= 1 ) {
                    for ( var k = 0; k < yRoots.length; k++ ) {
                        if ( Math.abs( xRoot - yRoots[k] ) < TOLERANCE ) {
                            var point = new Point(
                                c23x*(s*s*s)+(c22x*(s*s)+(c21x*(s)+(c20x))),
                                c23y*(s*s*s)+(c22y*(s*s)+(c21y*(s)+(c20y)))
                            );
                            locations.push( new CurveLocation( curve1, null, point, curve2 ) );
                            break checkRoots;
                        }
                    }
                }
            }
        }
    }
};


function findRoots( _p ){
    var n = _p.length-1, i, j;
    var z = 1, fz, dz;
    var roots = [], p = [], dp = [];
    var rootCount = 0, lastSign ,sign, pi, pj;
    lastSign = _p[0] > 0;
    // dp is the derivative of polynomial p w.r.t t
    for (i = 0, j = n; i <= n ; i++, j--){
        pi = _p[i];
        pj = _p[j];
        sign = pj > 0;
        if( i && pj && lastSign !== sign ){
            ++rootCount;
        }
        if( pj ) lastSign = sign;
        p.push( pi );
        dp.push( i * pi );
    }
    // According to Descartes' rule of signs, the number of real positive roots
    // is equal to or less than the change in signs of consecutive non-zero coeff.
    if( !rootCount )
        return roots;
    dp.shift();
    while( n-- ){
        for (i = 0; i <= MAX_ITERATE; i++) {
            fz = evaluateHorner( p, z );
            dfz = evaluateHorner( dp, z );
            dz = fz / dfz;
            // Check if we are done
            if (Math.abs(dz) < EPSILON){
                if( z >= 0 && z <= 1  ){
                    roots.push( z );
                }
                break;
            }
            // find the next approximation
            z = z - dz;
        }
        // z += EPSILON;
        // DEBUG: Should we exit if the root fails to converge?
        if( i > MAX_ITERATE ) break;
        // Divide p itself instead of creating a new array every time. Much faster this way!
        polyLongDivideSelf( p, z );
        dp.pop();
        for (i = 1; i <= n; i++)
            dp[i-1] = i * p[i];
    }
    return roots;
}

// TODO: Find a home for this one
/**
 * This method is analogous to paperjs#PathItem.getIntersections
 */
function getIntersections3( path1, path2 ){
    // First check the bounds of the two paths. If they don't intersect,
    // we don't need to iterate through their curves.
    if (!path1.getBounds().touches(path2.getBounds()))
        return [];
    var locations = [],
        curves1 = path1.getCurves(),
        curves2 = path2.getCurves(),
        length2 = curves2.length,
        values2 = [], ix, roots, len, implic;
    for (var i = 0; i < length2; i++)
        values2[i] = curves2[i].getValues();
    for (var i = 0, l = curves1.length; i < l; i++) {
        var curve1 = curves1[i],
            values1 = curve1.getValues();
        var v1Linear = Curve.isLinear(values1);
        // if(!v1Linear){
        //     implic = bernsteinToImplicit( values1 );
        // }
        for (var j = 0; j < length2; j++){
            value2 = values2[j];
            var v2Linear = Curve.isLinear(value2);
            if( v1Linear && v2Linear ){
                _getLineLineIntersection(values1, value2, curve1, curves2[j], locations);
            } else if ( v1Linear || v2Linear ){
                _getCurveLineIntersection(values1, value2, curve1, curves2[j], locations);
            } else {
                intersectCurveCurve(values1, value2, curve1, curves2[j], locations);
                // ix = getIntersectEquation( implic, value2 );
                // roots = findRoots( ix );
                // // str = "";
                // // len = 9;
                // // while( len-- )
                // //     str += "+ t^"+(len+1)+" "+ix[len]+" ";
                // len = roots.length;
                // // The roots are based on curve2, since we used curve1's implicit form
                // while( len-- ){
                //     var loc = new CurveLocation( curve1, null, curves2[j].getPointAt(roots[len], true), curves2[j] );
                //     if( loc.getParameter() )
                //         locations.push( loc );
                // }
                // console.log(str)
            }
        }
    }
    return locations;
}

