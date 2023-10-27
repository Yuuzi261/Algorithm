#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>

#define PI acos(-1);
#define Earth_R 6371.0;     // Radius of the Earth (km)

using namespace std;

class Coordinate {
    public:
        Coordinate();
        Coordinate(double, double);
        void set_coordinate(double, double);
        double get_longitude();
        double get_latitude();
    private:
        double longitude, latitude;
};

class POI {
    public:
        POI();
        POI(double, double);
        POI(double, double, string);
        void set_coordinate(double, double);
        void set_coordinate(Coordinate);
        double get_longitude();
        double get_latitude();
        string get_name();
        Coordinate get_coordinate();
        Coordinate get_cartesian_coordinate();
    private:
        Coordinate coordinate, cartesian_coordinate;
        string name;
        void cal_cartesian_coordinate();
};

class Edge {
    public:
        Edge(POI, POI);
        void set_edge(POI, POI);
        vector<POI> get_edge();
        double get_distance();
        double get_equ_c();
        double equ_substitute(double, double);
        double equ_substitute(POI);
    private:
        POI poi1, poi2;
        double distance;
        double equ_slope, equ_const, equ_longitude_cof, equ_latitude_cof, equ_c;
        void cal_equ();
};

class Map {
    public:
        Map();
        void add_POI(double, double);
        void add_POI(double, double, string);
        void add_Convex_Hull_edge(Edge);
        double get_POI_longitude(int);
        double get_POI_latitude(int);
        string get_POI_name(int);
        Coordinate get_POI_coordinate(int);
        POI get_POI(int);
        Edge get_CH_edge(int);
        int get_POIs_length();
        int get_CH_edges_length();
    private:
        vector<POI> POIs;
        vector<Edge> Convex_Hull_edges;
        int POIs_length, CH_edges_length;
};

// algorithms
void BruteForceClostestPoints(Map&, Coordinate&, Coordinate&);
void BruteForceConvexHull(Map&);

// other functions
double cal_Convex_Hull_area(Map&);

// tools
double to_radians(double);
double cal_distance(POI, POI);
double cal_Cartesian_distance(POI, POI);

int main(void) {

    Map map;

    ifstream ifs;

    int number_POIs;

    ifs.open("td.txt");
    ifs >> number_POIs;

    for(int i = 0;i < number_POIs;i++) {
        double longitude, latitude;
        string name;
        ifs >> longitude >> latitude >> name;
        map.add_POI(longitude, latitude, name);
    }

    ifs.close();

    cout << "[POIs List]" << endl;
    for(int i = 0;i < number_POIs;i++) {
        cout << map.get_POI_name(i) << ": (" << map.get_POI_longitude(i) << ", " << map.get_POI_latitude(i) << ')' << endl;
    }
    cout << endl;

    Coordinate coo1, coo2;
    BruteForceClostestPoints(map, coo1, coo2);
    cout << "Closest Pair: (" << coo1.get_longitude() << ", " << coo1.get_latitude() << "), (" << coo2.get_longitude() << ", " << coo2.get_latitude() << ')' << endl;

    BruteForceConvexHull(map);
    
    cout << endl << "[Convex Hull Edge List]" << endl;
    double total_distance = 0;
    for(int i = 0;i < map.get_CH_edges_length();i++) {
        Edge edge = map.get_CH_edge(i);
        vector<POI>pois = edge.get_edge();
        Coordinate coo1 = pois[0].get_coordinate(), coo2 = pois[1].get_coordinate();
        cout << pois[0].get_name() << ": (" << coo1.get_longitude() << ", " << coo1.get_latitude() << ") --- " << pois[1].get_name() << ": (" << coo2.get_longitude() << ", " << coo2.get_latitude() << ')' << endl;
        total_distance += edge.get_distance();
    }
    cout << endl << "Total Distance: " << total_distance << " km" << endl;
    cout << "Total Area: " << cal_Convex_Hull_area(map) << " km^2" << endl;

    return 0;

}

// Coordinate class
Coordinate::Coordinate() { this->longitude = 0; this->latitude = 0; }
Coordinate::Coordinate(double longitude, double latitude) { set_coordinate(longitude, latitude); }
void Coordinate::set_coordinate(double longitude, double latitude) { this->longitude = longitude; this->latitude = latitude; }
double Coordinate::get_longitude() { return this->longitude; }
double Coordinate::get_latitude() { return this->latitude; }

// POI class
POI::POI() { this->coordinate = Coordinate(); this->cartesian_coordinate = Coordinate(); this->name = ""; }
POI::POI(double longitude, double latitude) { this->coordinate = Coordinate(longitude, latitude); this->cal_cartesian_coordinate(); this->name = ""; }
POI::POI(double longitude, double latitude, string name) { this->coordinate = Coordinate(longitude, latitude); this->cal_cartesian_coordinate(); this->name = name; }
void POI::set_coordinate(double longitude, double latitude) { this->set_coordinate(longitude, latitude); }
void POI::set_coordinate(Coordinate coordinate) { this->coordinate = coordinate; }
double POI::get_longitude() { return this->coordinate.get_longitude(); }
double POI::get_latitude() { return this->coordinate.get_latitude(); }
string POI::get_name() { return this->name; }
Coordinate POI::get_coordinate() { return this->coordinate; }
Coordinate POI::get_cartesian_coordinate() { return this->cartesian_coordinate; }
void POI::cal_cartesian_coordinate() {
    double longitude = this->get_longitude(), latitude = this->get_latitude();
    double x = cos(to_radians(latitude)) * cos(to_radians(longitude)) * Earth_R;
    double y = cos(to_radians(latitude)) * sin(to_radians(longitude)) * Earth_R;
    this->cartesian_coordinate = Coordinate(x, y);
}

// Edge class
Edge::Edge(POI poi1, POI poi2) { set_edge(poi1, poi2); }
void Edge::set_edge(POI poi1, POI poi2) { this->poi1 = poi1; this->poi2 = poi2; this->distance = cal_Cartesian_distance(poi1, poi2); this->cal_equ(); }
vector<POI> Edge::get_edge() { vector<POI> edge = {this->poi1, this->poi2}; return edge; }
double Edge::get_distance() { return this->distance; }
double Edge::get_equ_c() { return this->equ_c; }
void Edge::cal_equ() {
    double delta_longitude = poi1.get_longitude() - poi2.get_longitude(), delta_latitude = poi1.get_latitude() - poi2.get_latitude();
    this->equ_slope = (delta_latitude) / (delta_longitude);
    this->equ_const = poi1.get_latitude() - this->equ_slope * poi1.get_longitude();
    this->equ_longitude_cof = -1 * delta_latitude;
    this->equ_latitude_cof = delta_longitude;
    this->equ_c = delta_longitude * this->equ_const;
}
double Edge::equ_substitute(double longitude, double latitude) { return this->equ_longitude_cof * longitude + this->equ_latitude_cof * latitude; }
double Edge::equ_substitute(POI poi) { return this->equ_longitude_cof * poi.get_longitude() + this->equ_latitude_cof * poi.get_latitude(); }

// Map class
Map::Map() { this->POIs_length = 0; this->CH_edges_length = 0; }
void Map::add_POI(double longitude, double latitude) { this->POIs.push_back(POI(longitude, latitude)); this->POIs_length++; }
void Map::add_POI(double longitude, double latitude, string name) { this->POIs.push_back(POI(longitude, latitude, name)); this->POIs_length++; }
void Map::add_Convex_Hull_edge(Edge edge) { this->Convex_Hull_edges.push_back(edge); this->CH_edges_length++; }
double Map::get_POI_longitude(int index) { return this->POIs[index].get_longitude(); }
double Map::get_POI_latitude(int index) { return this->POIs[index].get_latitude(); }
string Map::get_POI_name(int index) { return this->POIs[index].get_name(); }
Coordinate Map::get_POI_coordinate(int index) { return this->POIs[index].get_coordinate(); }
POI Map::get_POI(int index) { return this->POIs[index]; }
Edge Map::get_CH_edge(int index) { return this->Convex_Hull_edges[index]; }
int Map::get_POIs_length() { return this->POIs_length; }
int Map::get_CH_edges_length() { return this->CH_edges_length; }

// algorithms
void BruteForceClostestPoints(Map &map, Coordinate &coo1, Coordinate &coo2) {

    double dmin = INT32_MAX;

    for(int i = 0;i < map.get_POIs_length() - 1;i++) {
        for(int j = i+1;j < map.get_POIs_length();j++) {
            double distance = cal_distance(map.get_POI(i), map.get_POI(j));
            if(distance < dmin) {
                dmin = distance;
                coo1 = map.get_POI_coordinate(i);
                coo2 = map.get_POI_coordinate(j);
            }
        }
    }

}

void BruteForceConvexHull(Map &map) {

    for(int i = 0;i < map.get_POIs_length() - 1;i++) {
        for(int j = i+1;j < map.get_POIs_length();j++) {
            Edge edge(map.get_POI(i), map.get_POI(j));
            double less = 0, more = 0;
            bool is_convex_hull_edge = true;
            for(int k = 0; k < map.get_POIs_length();k++) {
                if(k != i && k != j) {
                    double c = edge.equ_substitute(map.get_POI(k));
                    if(c > edge.get_equ_c()) more++;
                    else if(c < edge.get_equ_c()) less++;
                    if(less != 0 && more != 0) {
                        is_convex_hull_edge = false;
                        break;
                    }
                }
            }
            if(is_convex_hull_edge) { map.add_Convex_Hull_edge(edge); break; }
        }
    }

    map.add_Convex_Hull_edge(Edge(map.get_CH_edge(map.get_CH_edges_length()-1).get_edge()[1], map.get_CH_edge(0).get_edge()[0]));

}

// other functions
double cal_Convex_Hull_area(Map &map) {

    double area = 0;

    for(int i = 0;i < map.get_CH_edges_length();i++) {
        vector<POI> pois = map.get_CH_edge(i).get_edge();
        Coordinate cacoo1 = pois[0].get_cartesian_coordinate(), cacoo2 = pois[1].get_cartesian_coordinate();
        area += cacoo1.get_longitude() * cacoo2.get_latitude() - cacoo2.get_longitude() * cacoo1.get_latitude();
    }

    return abs(area) / 2.0;
}

// tools
double to_radians(double degrees) { return degrees / 180.0 * PI; }
double cal_distance(POI poi1, POI poi2) { return sqrt(pow(poi1.get_longitude() - poi2.get_longitude(), 2) + pow(poi1.get_latitude() - poi2.get_latitude(), 2)); }
double cal_Cartesian_distance(POI poi1, POI poi2) {
    Coordinate cacoo1 = poi1.get_cartesian_coordinate(), cacoo2 = poi2.get_cartesian_coordinate();
    return sqrt(pow(cacoo1.get_longitude() - cacoo2.get_longitude(), 2) + pow(cacoo1.get_latitude() - cacoo2.get_latitude(), 2));
}

