
var EPSILON = 10e-12;
var TOLERANCE = 10e-6;

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
function polyLongDivideHorner( v, d ){
    var n = v.length-1, h = v[n], dh, i;
    var v1 = [];
    for( i = n-1; i >= 0; i-- ){
        dh = d * h;
        v1.push( dh );
        h= dh + v[i];
    }
    return v1;
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

function findRoot( v ){
    var n = v.length-1, i;
}
