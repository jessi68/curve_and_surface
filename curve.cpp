#include "curve.h"
#include "extra.h"
#ifdef WIN32
#include <windows.h>
#endif
#include <GL/gl.h>
using namespace std;

namespace
{
    // Approximately equal to.  We don't want to use == because of
    // precision issues with floating point.
    inline bool approx( const Vector3f& lhs, const Vector3f& rhs )
    {
        const float eps = 1e-8f;
        return ( lhs - rhs ).absSquared() < eps;
    }
    int DP_SIZE = 1001;
    int dpoints[1001][1001];
}

void initialize()
{
    for (int i = 0; i < DP_SIZE; i++) {
        for (int j = 0; j < DP_SIZE; j++) {
            dpoints[i][j] = -1;
        }
    }
}

int binomialCoefficient(int n, int k) {
    if (dpoints[n][k] != -1) {
        return dpoints[n][k];
    }

    if (k == 0 || k == n) {
        dpoints[n][k] = 1;
        return dpoints[n][k];
    }

    dpoints[n][k] = binomialCoefficient(n - 1, k - 1) + binomialCoefficient(n - 1, k);
    return dpoints[n][k];
}
    
Vector3f linearCombination(const vector< Vector3f >& points, float t) {
    Vector3f result = Vector3f(0, 0, 0);
    int degree = points.size() - 1;
    int currentDegree = degree;
    int currentOrder = 0;
    float coefficient;

    for (const auto& point : points) {
        coefficient = pow(1 - t, currentDegree)  * pow(t, currentOrder) * binomialCoefficient(degree, currentOrder);
        result += point * coefficient;
        currentDegree -= 1;
        currentOrder += 1;
    }

    return result;
}

Vector3f tangentOfBeizer(const vector< Vector3f >& points, float t) {
    Vector3f result = Vector3f(0, 0, 0);
    float coefficient;
    int size = points.size();
    int degree = size - 1;

    for (int i = 0; i <= degree - 1; i++) {
        result += binomialCoefficient(degree - 1, i) * size * (points[i + 1] - points[i]);
    }
    return result;
}

Vector3f normalOfBeizer(const vector< Vector3f >& points, float t) {
    Vector3f result = Vector3f(0, 0, 0);
    return result;
}

Curve evalBezier(const vector< Vector3f >& points, unsigned steps)
{
    // Check
    // 4, 7, 10, 13 
    //Curve Bezier(steps + 1);
    int orderOfCurve = points.size();
    cout << orderOfCurve << endl;
    cout << "is bezier called" << endl;
    if (orderOfCurve < 4 || orderOfCurve % 3 != 1)
    {
        cerr << "evalBezier must be called with 3n+1 control points." << endl;
        exit(0);
    }

    int seg;			
    //Determine number of segments if 4 or 
    //more control points
    if (points.size() == 4) {
        seg = 1;
    }
    else {
        seg = ((points.size() - 1) / 3);
    }

    Curve R((seg) * (steps)+1);

    Vector3f initialBinormal;		//Initial Binormal (arbitrary only at beginning)
    Vector3f recentBinormal;		//Most recent Binormal

    recentBinormal = Vector3f(0, 0, 1);
    int segment = 0;

    for (unsigned i = 0; i < orderOfCurve - 3; i += 3)
    {

        Matrix4f vertexMatrix;
        Matrix3f tangentMatrix;
        Matrix2f normalMatrix;
        Vector4f time;
        Matrix4f controlPoints;
        Matrix3f diffOfControlPoints;

        Vector4f vertexResult;
        int currentIndex;
        
        for (unsigned delta = 0; delta <= steps; delta += 1)
        {
            float t = float(delta) / steps;
            vertexMatrix = Matrix4f(1, -3, 3, -1,
                0, 3, -6, 3,
                0, 0, 3, -3,
                0, 0, 0, 1);
            tangentMatrix = Matrix3f(3, -6, 3,
                0, 6, -6,
                0, 0, 3);
            time = Vector4f(1, t, t * t, t * t * t);
            controlPoints = Matrix4f(points[i + 0][0], points[i + 1][0], points[i + 2][0], points[i + 3][0],
                points[i + 0][1], points[i + 1][1], points[i + 2][1], points[i + 3][1],
                points[i + 0][2], points[i + 1][2], points[i + 2][2], points[i + 3][2],
                0, 0, 0, 0);
            
            diffOfControlPoints = Matrix3f(points[i + 1] - points[i], points[i + 2] - points[i + 1], points[i + 3] - points[i + 2], true);
            // T is transpose Matrix 
            // in opengl column first so we have to take transpose matix 
            // (time * vertexMatrix * controlPoints)T = (controlPoints)T * (vertexMatrix)T * (time)T 
            vertexResult = (controlPoints * vertexMatrix * time);

            cout << i << " i " << endl;
            cout << steps << endl;
            cout << delta << endl;
            currentIndex = steps * segment + delta;
        
            if (currentIndex >= R.size()) {
                break;
            }
            R[currentIndex].V = vertexResult.xyz();
            R[currentIndex].T = (diffOfControlPoints * tangentMatrix * time.xyz()).normalized();
            // infection point 방지 위해 
            // 또한 직선운동일 때, N 은 0 이 되버림 
            R[currentIndex].N = Vector3f::cross(recentBinormal, R[currentIndex].T).normalized();
            R[currentIndex].B = Vector3f::cross(R[currentIndex].T, R[currentIndex].N).normalized();
            recentBinormal = R[currentIndex].B;
        }
        segment += 1;
    }

    return R;

}

Curve evalBspline( const vector< Vector3f >& P, unsigned steps )
{
    // Check
    if( P.size() < 4 )
    {
        cerr << "evalBspline must be called with 4 or more control points." << endl;
        exit( 0 );
    }

    // TODO:
    // It is suggested that you implement this function by changing
    // basis from B-spline to Bezier.  That way, you can just call
    // your evalBezier function.

    cerr << "\t>>> evalBSpline has been called with the following input:" << endl;

    cerr << "\t>>> Control points (type vector< Vector3f >): "<< endl;
    for( unsigned i = 0; i < P.size(); ++i )
    {
        cerr << "\t>>> " << P[i] << endl;
    }
    
    cerr << "\t>>> Steps (type steps): " << steps << endl;
    cerr << "\t>>> Returning empty curve." << endl;

    // Return an empty curve right now.
    return Curve();
}

Curve evalCircle( float radius, unsigned steps )
{
    // This is a sample function on how to properly initialize a Curve
    // (which is a vector< CurvePoint >).
    
    // Preallocate a curve with steps+1 CurvePoints
    Curve R( steps+1 );

    // Fill it in counterclockwise
    for( unsigned i = 0; i <= steps; ++i )
    {
        // step from 0 to 2pi
        float t = 2.0f * M_PI * float( i ) / steps;

        // Initialize position
        // We're pivoting counterclockwise around the y-axis
        R[i].V = radius * Vector3f( cos(t), sin(t), 0 );
        
        // Tangent vector is first derivative
        R[i].T = Vector3f( -sin(t), cos(t), 0 );
        
        // Normal vector is second derivative
        R[i].N = Vector3f( -cos(t), -sin(t), 0 );

        // Finally, binormal is facing up.
        R[i].B = Vector3f( 0, 0, 1 );
    }

    return R;
}

void drawCurve( const Curve& curve, float framesize )
{
    // Save current state of OpenGL
    glPushAttrib( GL_ALL_ATTRIB_BITS );

    // Setup for line drawing
    glDisable( GL_LIGHTING ); 
    glColor4f( 1, 1, 1, 1 );
    glLineWidth( 1 );
    
    // 


    
    cout << curve.size();
    cout << "draw curve" << endl;

    
    glBegin(GL_LINE_STRIP);
    for (unsigned i = 0; i < curve.size(); ++i)
    {
        cout << "curve coordinate";
        cout << curve[i].V.x() << " " << curve[i].V.y() << curve[i].V.z() << endl;
        glVertex(curve[i].V);
    }
    glEnd();
    

    // Draw coordinate frames if framesize nonzero
    cout << framesize << endl;
    cout << "currve" << endl;
    if( framesize != 0.0f )
    {
        Matrix4f M;

        for( unsigned i = 0; i < curve.size(); ++i )
        {
            M.setCol( 0, Vector4f( curve[i].N, 0 ) );
            M.setCol( 1, Vector4f( curve[i].B, 0 ) );
            M.setCol( 2, Vector4f( curve[i].T, 0 ) );
            M.setCol( 3, Vector4f( curve[i].V, 1 ) );
            
            glPushMatrix();
            glMultMatrixf( M );
            glScaled( framesize, framesize, framesize );
            glBegin( GL_LINES );
            glColor3f( 1, 0, 0 ); glVertex3d( 0, 0, 0 ); glVertex3d( 1, 0, 0 );
            glColor3f( 0, 1, 0 ); glVertex3d( 0, 0, 0 ); glVertex3d( 0, 1, 0 );
            glColor3f( 0, 0, 1 ); glVertex3d( 0, 0, 0 ); glVertex3d( 0, 0, 1 );
            glEnd();
            glPopMatrix();
        }
    }
    
    // Pop state
    glPopAttrib();
}

