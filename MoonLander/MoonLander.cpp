#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>

const double PI = 3.1415927;
const double DEGREES_TO_RAD = PI / 180.0;
const double RAD_TO_DEGREES = 180.0 / PI;

const double SHIP_RADIUS = 100.0;

class Point
{
public:
    Point()
        : x(0.0)
        , y(0.0)
        , z(0.0)
    {}

    Point(double _x, double _y, double _z = 0.0)
        : x(_x)
        , y(_y)
        , z(_z)
    {}

    double x;
    double y;
    double z;
};
 
class Line
{
public:
    Line()
        : pt()
    {}

    Line(const Point &p1, const Point &p2)
        : pt()
    {
        pt[0] = p1;
        pt[1] = p2;
    }

    bool IsFlat() const
    {
        return pt[0].y == pt[1].y;
    }

    double Slope() const
    {
        if ((pt[1].x - pt[0].x) != 0.0)
        {
            return (pt[1].y - pt[0].y) / (pt[1].x - pt[0].x);
        }
        else
        {
            return 0.0;
        }
    }

    Point pt[2];
};

 
double dot(Point c1, Point c2)
{
    return (c1.x * c2.x + c1.y * c2.y + c1.z * c2.z);
}
 
double norm(Point c1)
{
    return std::sqrt(dot(c1, c1));
}
 
double getShortestDistance(Line line1, Line line2)
{
    double EPS = 0.00000001;
 
    Point delta21;
    delta21.x = line1.pt[1].x - line1.pt[0].x;
    delta21.y = line1.pt[1].y - line1.pt[0].y;
    delta21.z = line1.pt[1].z - line1.pt[0].z;
 
    Point delta41;
    delta41.x = line2.pt[1].x - line2.pt[0].x;
    delta41.y = line2.pt[1].y - line2.pt[0].y;
    delta41.z = line2.pt[1].z - line2.pt[0].z;
 
    Point delta13;
    delta13.x = line1.pt[0].x - line2.pt[0].x;
    delta13.y = line1.pt[0].y - line2.pt[0].y;
    delta13.z = line1.pt[0].z - line2.pt[0].z;
 
    double a = dot(delta21, delta21);
    double b = dot(delta21, delta41);
    double c = dot(delta41, delta41);
    double d = dot(delta21, delta13);
    double e = dot(delta41, delta13);
    double D = a * c - b * b;
 
    double sc, sN, sD = D;
    double tc, tN, tD = D;
 
    if (D < EPS)
    {
        sN = 0.0;
        sD = 1.0;
        tN = e;
        tD = c;
    }
    else
    {
        sN = (b * e - c * d);
        tN = (a * e - b * d);
        if (sN < 0.0)
        {
            sN = 0.0;
            tN = e;
            tD = c;
        }
        else if (sN > sD)
        {
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }
 
    if (tN < 0.0)
    {
        tN = 0.0;
 
        if (-d < 0.0)
            sN = 0.0;
        else if (-d > a)
            sN = sD;
        else
        {
            sN = -d;
            sD = a;
        }
    }
    else if (tN > tD)
    {
        tN = tD;
        if ((-d + b) < 0.0)
            sN = 0;
        else if ((-d + b) > a)
            sN = sD;
        else
        {
            sN = (-d + b);
            sD = a;
        }
    }
 
    if (std::abs(sN) < EPS) sc = 0.0;
    else sc = sN / sD;
    if (std::abs(tN) < EPS) tc = 0.0;
    else tc = tN / tD;
 
    Point dP;
    dP.x = delta13.x + (sc * delta21.x) - (tc * delta41.x);
    dP.y = delta13.y + (sc * delta21.y) - (tc * delta41.y);
    dP.z = delta13.z + (sc * delta21.z) - (tc * delta41.z);
 
    return std::sqrt(dot(dP, dP));
}

class Vector
{
public:
    double angle;
    double mag;

    friend Vector operator+(const Vector &v1, const Vector &v2)
    {
        if (v1.mag == 0)
        {
            return v2;
        }
        else if (v2.mag == 0)
        {
            return v1;
        }
        else
        {
            double v1_x = v1.mag * std::cos(v1.angle * DEGREES_TO_RAD);
            double v1_y = v1.mag * std::sin(v1.angle * DEGREES_TO_RAD);
            double v2_x = v2.mag * std::cos(v2.angle * DEGREES_TO_RAD);
            double v2_y = v2.mag * std::sin(v2.angle * DEGREES_TO_RAD);
            double v_x = v1_x + v2_x;
            double v_y = v1_y + v2_y;

            Vector v;
            v.angle = std::atan(v_y / v_x) * RAD_TO_DEGREES;
            v.mag = std::sqrt(v_x * v_x + v_y * v_y);

            return v;
        }
    }
};

class PIDController
{
public:
    PIDController()
    {}

private:
    double m_kp;
    double m_ki;
    double m_kd;

    double m_error;
    double m_error_prev;
    double m_target;
};

class Autopilot
{
public:
    Autopilot()
        : m_in()
        , m_surface_points()
        , m_surface_lines()
        , m_ship_loc()
        , m_flight_path()
        , m_is_initialized(false)
        , m_dmd_angle(0)
        , m_dmd_thrust(0)
        , m_target_alt(0.0)
        , m_debug_mode(false)
    {}

    struct Input
    {
        int x;
        int y;
        int h_speed;   // the horizontal speed (in m/s), can be negative.
        int v_speed;   // the vertical speed (in m/s), can be negative.
        int fuel;      // the quantity of remaining fuel in liters.
        int rotate;    // the rotation angle in degrees (-90 to 90).
        int power;     // the thrust power (0 to 4).
    };

    void SetDebugMode(bool val) { m_debug_mode = val; }

    Input & In() { return m_in; }

    void AddSurfacePoint(const Point &p)
    {
        m_surface_points.push_back(p);
    }
    
    void Update()
    {
        m_ship_loc = Point(m_in.x, m_in.y);
        if (!m_is_initialized)
        {
            Initialize();
        }

        
        if(CheckLandingMode())
        {
            PerformLanding();
        }
        else
        {
            ComputeAltitudeTarget();
            ComputeHorizontalSpeedTarget();
            ComputeThrustCommands();
        }

        std::cout << m_dmd_angle << " " << m_dmd_thrust << std::endl;
    }
private:
    bool CheckLandingMode()
    {
        return m_landing_surface.pt[0].x < m_ship_loc.x && m_ship_loc.x < m_landing_surface.pt[1].x && std::abs(m_in.h_speed) < 15.0;
    }

    void ComputeThrustCommands()
    {
        double alt_error = m_target_alt - m_ship_loc.y;
        double h_speed_error = m_target_h_speed - m_in.h_speed;
        double v_speed_error = -35.0 - m_in.v_speed;

        if (m_debug_mode)
        {
            std::cerr << "Alt Error: " << alt_error << std::endl;
            std::cerr << "Vspeed Error: " << v_speed_error << std::endl;
            std::cerr << "Speed Error: " << h_speed_error << std::endl;
        }

        Vector alt_vector;
        alt_vector.angle = 0.0;
        alt_vector.mag = std::max(0.2 * alt_error, 0.0);
        std::cerr << "Alt Vector : " << alt_vector.angle << "°, " << alt_vector.mag << std::endl;
        
        

        Vector v_speed_vector;
        v_speed_vector.angle = 0.0;
        v_speed_vector.mag = std::max(2.0 * v_speed_error, 0.0);
        std::cerr << "Vspeed Vector : " << v_speed_vector.angle << "°, " << v_speed_vector.mag << std::endl;


        Vector h_speed_vector;
        h_speed_vector.mag = 0.25 * h_speed_error;
        h_speed_vector.angle = (h_speed_vector.mag > 0.0) ? -70.0 : 70.0;
        h_speed_vector.mag = std::abs(h_speed_vector.mag);
        std::cerr << "Speed Vector : " << h_speed_vector.angle << "°, " << h_speed_vector.mag << std::endl;

        Vector thrust_vect = alt_vector + h_speed_vector + v_speed_vector;
        
        std::cerr << "Thrust Vector : " << thrust_vect.angle << "°, " << thrust_vect.mag << std::endl;

        m_dmd_angle =  static_cast<int>(thrust_vect.angle);
        if(std::abs(m_in.rotate - m_dmd_angle) < 5.0)
        {
            m_dmd_thrust = static_cast<int>(std::max(std::min(thrust_vect.mag, 4.0), 0.0));
        }
        else
        {
            m_dmd_thrust = 0;
        }
    }

    void ComputeAltitudeTarget()
    {
        int curr_path_tgt_idx = 0;
        for (auto &pt : m_flight_path)
        {
            if (pt.x < m_ship_loc.x)
            {
                ++curr_path_tgt_idx;
            }
            else 
            {
                break;
            }
        }
        Line curr_flight_path;
        if (curr_path_tgt_idx > 0)
        {
            curr_flight_path.pt[0] = m_flight_path[curr_path_tgt_idx - 1];
            curr_flight_path.pt[1] = m_flight_path[curr_path_tgt_idx];
        }
        else
        {
            curr_flight_path.pt[0] = m_flight_path[0];
            curr_flight_path.pt[1] = m_flight_path[1];
        }

        double m = curr_flight_path.Slope();
        double b = curr_flight_path.pt[0].y - (m * curr_flight_path.pt[0].x);
        m_target_alt = m * m_ship_loc.x + b;

        if (m_debug_mode)
        {
            std::cerr << "Target Altitude: " << m_target_alt << std::endl;
        }
    }

    void ComputeHorizontalSpeedTarget()
    {
        double lz_x_diff = m_flight_path.back().x - m_ship_loc.x;
        
        if (lz_x_diff < 0.0)
        {
            m_target_h_speed = std::min(lz_x_diff * 0.01, static_cast<double>(m_initial_h_speed));
        }
        else
        {
            m_target_h_speed = std::max(lz_x_diff * 0.01, static_cast<double>(m_initial_h_speed));
        }
        
        if (m_debug_mode)
        {
            std::cerr << "LZ X Diff: " << lz_x_diff << std::endl;
            std::cerr << "Speed Target: " << m_target_h_speed << std::endl;
        }
    }

    void Initialize()
    {
        m_initial_h_speed = m_in.h_speed;
        for (int i = 0; i < m_surface_points.size() - 1; ++i)
        {
            m_surface_lines.push_back(Line(m_surface_points[i], m_surface_points[i + 1]));
        }

        if (m_debug_mode)
        {
            for (auto & pt : m_surface_points)
            {
                //std::cerr << "Surface Pt: " << pt.x << ", " << pt.y << std::endl;
            }
        }

        FindFlightPath();

        m_is_initialized = true;
    }
    
    void PerformLanding()
    {
        std::cerr << "Landing Mode..." << std::endl;
        double v_speed_error = -35.0 - m_in.v_speed;

        if (m_debug_mode)
        {
            std::cerr << "Vspeed Error: " << v_speed_error << std::endl;
        }

        Vector alt_vector;
        alt_vector.angle = 0.0;
        alt_vector.mag = std::max(2.0 * v_speed_error, 0.0);
        std::cerr << "Vspeed Vector : " << alt_vector.angle << "°, " << alt_vector.mag << std::endl;

        

        Vector thrust_vect = alt_vector;
        
        std::cerr << "Thrust Vector : " << thrust_vect.angle << "°, " << thrust_vect.mag << std::endl;

        m_dmd_angle =  static_cast<int>(thrust_vect.angle);
        m_dmd_thrust = static_cast<int>(std::max(std::min(thrust_vect.mag, 4.0), 0.0));
    }

    void FindFlightPath()
    {
        // Find landing target;
        Point lz;
        Point flare_pt;
        for (auto & ln : m_surface_lines)
        {
            if (ln.IsFlat())
            {
                lz.x = ln.pt[0].x + (ln.pt[1].x - ln.pt[0].x) / 2.0;
                lz.y = ln.pt[0].y;
                flare_pt.x = ln.pt[0].x + 600.0;
                flare_pt.y = ln.pt[1].y + 100.0;
                m_landing_surface = ln;
            }
        }

        m_flight_path.push_back(m_ship_loc);

        Line flight_path_segment(m_ship_loc, flare_pt);

        for (auto &ln : m_surface_lines)
        {
            if (getShortestDistance(flight_path_segment, ln) <= SHIP_RADIUS && !ln.IsFlat())
            {
                flight_path_segment.pt[1].x = ln.pt[1].x;
                flight_path_segment.pt[1].y = ln.pt[1].y + 3.0 * SHIP_RADIUS;
                double slope = flight_path_segment.Slope();
                if (std::abs(slope) > 0.5)
                {
                    flight_path_segment.pt[1].y = 0.2 * (ln.pt[1].x - ln.pt[0].x) + ln.pt[0].y;
                }
                m_flight_path.push_back(flight_path_segment.pt[1]);
                flight_path_segment = Line(m_flight_path.back(), flare_pt);
            }
        }

        // We didnt find any additional waypoints...
        if (m_flight_path.size() == 1)
        {
            m_flight_path.push_back(flight_path_segment.pt[1]);
        }
        else
        {
            m_flight_path.push_back(flare_pt);
        }

        for (auto & ln : m_flight_path)
        {
            std::cerr << "Flight Path: " << ln.x << ", " << ln.y << std::endl;
        }
        
        m_lz_loc = lz;

    }
    Input m_in;
    std::vector<Point> m_surface_points;
    std::vector<Line> m_surface_lines;
    Line m_landing_surface;
    Point m_ship_loc;
    Point m_lz_loc;
    std::vector<Point> m_flight_path;
    double m_target_alt;
    double m_target_h_speed;
    int m_dmd_angle;
    int m_dmd_thrust;
    bool m_is_initialized;
    bool m_debug_mode;
    int m_initial_h_speed;
};

void test()
{
    Autopilot ap;
    ap.SetDebugMode(true);
    
    ap.AddSurfacePoint(Point(   0, 1000));
    ap.AddSurfacePoint(Point( 300, 1500));
    ap.AddSurfacePoint(Point( 350, 1400));
    ap.AddSurfacePoint(Point( 500, 2000));
    ap.AddSurfacePoint(Point( 800, 1800));
    ap.AddSurfacePoint(Point(1000, 2500));
    ap.AddSurfacePoint(Point(1200, 2100));
    ap.AddSurfacePoint(Point(1500, 2400));
    ap.AddSurfacePoint(Point(2000, 1000));
    ap.AddSurfacePoint(Point(2200,  500));
    ap.AddSurfacePoint(Point(2500,  100));
    ap.AddSurfacePoint(Point(2900,  800));
    ap.AddSurfacePoint(Point(3000,  500));
    ap.AddSurfacePoint(Point(3200, 1000));
    ap.AddSurfacePoint(Point(3500, 2000));
    ap.AddSurfacePoint(Point(3800,  800));
    ap.AddSurfacePoint(Point(4000,  200));
    ap.AddSurfacePoint(Point(5000,  200));
    ap.AddSurfacePoint(Point(5500, 1500));
    ap.AddSurfacePoint(Point(6999, 2800));

    ap.In().x = 500;
    ap.In().y = 2700;
    ap.Update();

    ap.In().x = 1100;
    ap.In().y = 2500;
    ap.Update();
}

void execute()
{
    Autopilot ap;
    ap.SetDebugMode(true);
    int surfaceN; // the number of points used to draw the surface of Mars.
    std::cin >> surfaceN; std::cin.ignore();
    for (int i = 0; i < surfaceN; i++) {
        int landX; // X coordinate of a surface point. (0 to 6999)
        int landY; // Y coordinate of a surface point. By linking all the points together in a sequential fashion, you form the surface of Mars.
        std::cin >> landX >> landY; std::cin.ignore();
        
        ap.AddSurfacePoint(Point(landX, landY));
    }

    // game loop
    while (1) {
        std::cin >> ap.In().x >> ap.In().y >> ap.In().h_speed >> ap.In().v_speed >> ap.In().fuel >> ap.In().rotate >> ap.In().power; std::cin.ignore();
        ap.Update();
    }
}

int main()
{
    //test();
    execute();
    return 0;
}

