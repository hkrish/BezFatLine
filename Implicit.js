
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
        v[i] = dh;
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
function bernseinToImplicit( v ){
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

    return [ A, B, C, D, E, F, G, H, I, J ];
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
     * get a 9 degree polynomial, roots can be calculated by evaluating
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

        return [ t0, t1, t2, t3, t4, t5, t6, t7, t8, t9 ];
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
    // dp -> d( p )/dt
    for (i = 0; i <= n; i++){
        p.push( _p[i] );
        dp.push( i * p[i] );
    }
    dp.shift();
    i = 0;
    while( n-- ){
        for (i = 0; i <= MAX_ITERATE; i++) {
            fz = evaluateHorner( p, z );
            dz = fz / evaluateHorner( dp, z );
            // Check if we are done
            if (Math.abs(dz) < EPSILON){
                if( z >= 0 && z <= 1 && evaluateHorner(_p, z) < TOLERANCE ){
                    roots.push( z );
                }
                break;
            }
            // find the next approximation
            z = z - dz;
        }
        // Divide p itself instead of creating a new array every time. Much faster this way!
        polyLongDivideSelf( p, z );
        dp.pop();
        for (i = 1; i <= n; i++)
            dp[i-1] = i * p[i];
    }
    return roots;
}
