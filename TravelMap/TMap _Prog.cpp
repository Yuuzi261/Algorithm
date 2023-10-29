#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <windows.h>

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
        void set_isCH_POI(bool);
        double get_longitude();
        double get_latitude();
        double get_x();
        double get_y();
        string get_name();
        bool get_isCH_POI();
        Coordinate get_coordinate();
        Coordinate get_cartesian_coordinate();
    private:
        Coordinate coordinate, cartesian_coordinate;
        string name;
        bool isCH_POI;
        void cal_cartesian_coordinate();
};

class Edge {
    public:
        Edge();
        Edge(POI, POI);
        void set_edge(POI, POI);
        vector<POI> get_edge();
        double get_distance();
        double get_equ_c();
        double equ_substitute(double, double);
        double equ_substitute(POI);
        double cal_distance_with_poi(POI);
    private:
        POI poi1, poi2;
        double distance;
        double equ_slope, equ_const, equ_x_cof, equ_y_cof, equ_c;
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
        vector<int> get_shortest_route();
        vector<POI> get_CH_shortest_route();
        double get_min_distance();
        int get_POIs_length();
        int get_CH_edges_length();
        // algorithms
        void BruteForceClostestPoints(POI&, POI&);
        void BruteForceConvexHull();
        void BruteForceTSP();
        void ConvexHullTSP();
    private:
        vector<POI> POIs;
        vector<Edge> Convex_Hull_edges;
        vector<int> shortest_route;
        vector<POI> CH_shortest_route;
        double min_distance;
};

// other functions
double cal_Convex_Hull_area(Map&);

// tools
double to_radians(double);
double cal_Cartesian_distance(POI, POI);
void sort_with_distance(vector<POI>&, POI);

int main(void) {

    LARGE_INTEGER freq, start, end;
    QueryPerformanceFrequency(&freq);

    Map map;

    ifstream ifs;

    int number_POIs;

    string map_filename;
    cout << "Enter the filename of Map: ";
    getline(cin, map_filename);

    ifs.open(map_filename + ".in");
    ifs >> number_POIs;

    for(int i = 0;i < number_POIs;i++) {
        double longitude, latitude;
        string name;
        ifs >> longitude >> latitude;
        getline(ifs, name);
        map.add_POI(longitude, latitude, string().assign(name, 1, name.size()));
    }

    ifs.close();

    cout << endl << "[POIs List]" << endl;
    for(int i = 0;i < number_POIs;i++) {
        cout << map.get_POI_name(i) << ": (" << map.get_POI_longitude(i) << ", " << map.get_POI_latitude(i) << ')' << endl;
    }
    cout << endl;

    POI poi1, poi2;
    map.BruteForceClostestPoints(poi1, poi2);
    cout << "[Closest Pair]" << endl << poi1.get_name() << ": (" << poi1.get_longitude() << ", " << poi1.get_latitude() << "), " << poi2.get_name() << ": (" << poi2.get_longitude() << ", " << poi2.get_latitude() << ')' << endl;
    cout << "Distance: " << cal_Cartesian_distance(poi1, poi2) << " km" << endl;

    map.BruteForceConvexHull();
    
    cout << endl << "[Convex Hull Edge List]" << endl;
    double total_distance = 0;
    for(int i = 0;i < map.get_CH_edges_length();i++) {
        Edge edge = map.get_CH_edge(i);
        vector<POI>pois = edge.get_edge();
        Coordinate coo1 = pois[0].get_coordinate(), coo2 = pois[1].get_coordinate();
        cout << pois[0].get_name() << ": (" << coo1.get_longitude() << ", " << coo1.get_latitude() << ") --> " << pois[1].get_name() << ": (" << coo2.get_longitude() << ", " << coo2.get_latitude() << ')' << endl;
        total_distance += edge.get_distance();
    }
    cout << endl << "Total Distance: " << total_distance << " km" << endl;
    cout << "Total Area: " << cal_Convex_Hull_area(map) << " km^2" << endl;

    QueryPerformanceCounter(&start);
    map.BruteForceTSP();
    QueryPerformanceCounter(&end);
    double BF_elapsed = (double)(end.QuadPart - start.QuadPart) / freq.QuadPart;

    cout << endl << "[TSP Problem]" << endl << "Brute Force TSP:" << " Elasped: " << BF_elapsed << " s" << endl;
    vector<int> shortest_route = map.get_shortest_route();
    for(int i = 0;i < shortest_route.size()-1;i++) {
        POI poi = map.get_POI(shortest_route[i]);
        cout << poi.get_name() << " --> ";
    }
    cout << map.get_POI(shortest_route[shortest_route.size()-1]).get_name() << endl;
    cout << "Total Distance: " << map.get_min_distance() << " km" << endl;

    QueryPerformanceCounter(&start);
    map.ConvexHullTSP();
    QueryPerformanceCounter(&end);
    double CH_elapsed = (double)(end.QuadPart - start.QuadPart) / freq.QuadPart * 1000;

    cout << endl << "Convex Hull TSP:" << " Elasped: " << CH_elapsed << " ms" << endl;
    vector<POI> CH_shortest_route = map.get_CH_shortest_route();
    double CH_min_distance = 0;
    for(int i = 0;i < CH_shortest_route.size()-1;i++){
        cout << CH_shortest_route[i].get_name() << " --> ";
        CH_min_distance += cal_Cartesian_distance(CH_shortest_route[i], CH_shortest_route[i+1]);
    }
    cout << CH_shortest_route[CH_shortest_route.size()-1].get_name() << endl << "Total Distance: " << CH_min_distance << " km" << endl << endl;

    cout << "Distance Ratio = " << CH_min_distance / map.get_min_distance() << endl;
    cout << "Execution Time Ratio = " << CH_elapsed / BF_elapsed << endl;

    return 0;

}

// Coordinate class
Coordinate::Coordinate() { this->longitude = 0; this->latitude = 0; }
Coordinate::Coordinate(double longitude, double latitude) { set_coordinate(longitude, latitude); }
void Coordinate::set_coordinate(double longitude, double latitude) { this->longitude = longitude; this->latitude = latitude; }
double Coordinate::get_longitude() { return this->longitude; }
double Coordinate::get_latitude() { return this->latitude; }

// POI class
POI::POI() { this->coordinate = Coordinate(); this->cartesian_coordinate = Coordinate(); this->name = ""; this->isCH_POI = false; }
POI::POI(double longitude, double latitude) { this->coordinate = Coordinate(longitude, latitude); this->cal_cartesian_coordinate(); this->name = ""; this->isCH_POI = false; }
POI::POI(double longitude, double latitude, string name) { this->coordinate = Coordinate(longitude, latitude); this->cal_cartesian_coordinate(); this->name = name; this->isCH_POI = false; }
void POI::set_coordinate(double longitude, double latitude) { this->set_coordinate(longitude, latitude); }
void POI::set_coordinate(Coordinate coordinate) { this->coordinate = coordinate; }
void POI::set_isCH_POI(bool isCH_POI) { this->isCH_POI = isCH_POI; }
double POI::get_longitude() { return this->coordinate.get_longitude(); }
double POI::get_latitude() { return this->coordinate.get_latitude(); }
double POI::get_x() { return this->cartesian_coordinate.get_longitude(); }
double POI::get_y() { return this->cartesian_coordinate.get_latitude(); }
string POI::get_name() { return this->name; }
bool POI::get_isCH_POI() { return this->isCH_POI; }
Coordinate POI::get_coordinate() { return this->coordinate; }
Coordinate POI::get_cartesian_coordinate() { return this->cartesian_coordinate; }
void POI::cal_cartesian_coordinate() {
    double longitude = this->get_longitude(), latitude = this->get_latitude();
    double x = cos(to_radians(latitude)) * cos(to_radians(longitude)) * Earth_R;
    double y = cos(to_radians(latitude)) * sin(to_radians(longitude)) * Earth_R;
    this->cartesian_coordinate = Coordinate(x, y);
}

// Edge class
Edge::Edge() {}
Edge::Edge(POI poi1, POI poi2) { set_edge(poi1, poi2); }
void Edge::set_edge(POI poi1, POI poi2) { this->poi1 = poi1; this->poi2 = poi2; this->distance = cal_Cartesian_distance(poi1, poi2); this->cal_equ(); }
vector<POI> Edge::get_edge() { vector<POI> edge = {this->poi1, this->poi2}; return edge; }
double Edge::get_distance() { return this->distance; }
double Edge::get_equ_c() { return this->equ_c; }
void Edge::cal_equ() {
    double delta_x = poi1.get_x() - poi2.get_x(), delta_y = poi1.get_y() - poi2.get_y();
    this->equ_slope = (delta_y) / (delta_x);
    this->equ_const = poi1.get_y() - this->equ_slope * poi1.get_x();
    this->equ_x_cof = -1 * delta_y;
    this->equ_y_cof = delta_x;
    this->equ_c = delta_x * this->equ_const;
}
double Edge::equ_substitute(double x, double y) { return this->equ_x_cof * x + this->equ_y_cof * y; }
double Edge::equ_substitute(POI poi) { return this->equ_x_cof * poi.get_x() + this->equ_y_cof * poi.get_y(); }
double Edge::cal_distance_with_poi(POI poi) { return abs(this->equ_x_cof * poi.get_x() + this->equ_y_cof * poi.get_y() - this->equ_c) / sqrt(pow(this->equ_x_cof, 2) + pow(this->equ_y_cof, 2)); }

// Map class
Map::Map() {}
void Map::add_POI(double longitude, double latitude) { this->POIs.push_back(POI(longitude, latitude)); }
void Map::add_POI(double longitude, double latitude, string name) { this->POIs.push_back(POI(longitude, latitude, name)); }
void Map::add_Convex_Hull_edge(Edge edge) { this->Convex_Hull_edges.push_back(edge); }
double Map::get_POI_longitude(int index) { return this->POIs[index].get_longitude(); }
double Map::get_POI_latitude(int index) { return this->POIs[index].get_latitude(); }
string Map::get_POI_name(int index) { return this->POIs[index].get_name(); }
Coordinate Map::get_POI_coordinate(int index) { return this->POIs[index].get_coordinate(); }
POI Map::get_POI(int index) { return this->POIs[index]; }
Edge Map::get_CH_edge(int index) { return this->Convex_Hull_edges[index]; }
vector<int> Map::get_shortest_route() { return this->shortest_route; }
vector<POI> Map::get_CH_shortest_route() { return this->CH_shortest_route; }
double Map::get_min_distance() { return this->min_distance; }
int Map::get_POIs_length() { return this->POIs.size(); }
int Map::get_CH_edges_length() { return this->Convex_Hull_edges.size(); }

// algorithms
void Map::BruteForceClostestPoints(POI &poi1, POI &poi2) {

    double dmin = INT32_MAX;

    for(int i = 0;i < this->get_POIs_length() - 1;i++) {
        for(int j = i+1;j < this->get_POIs_length();j++) {
            double distance = cal_Cartesian_distance(this->POIs[i], this->POIs[j]);
            if(distance < dmin) {
                dmin = distance;
                poi1 = this->POIs[i];
                poi2 = this->POIs[j];
            }
        }
    }

}

void Map::BruteForceConvexHull() {

    for(int i = 0;i < this->get_POIs_length() - 1;i++) {
        for(int j = i+1;j < this->get_POIs_length();j++) {
            Edge edge(this->POIs[i], this->POIs[j]);
            double less = 0, more = 0;
            bool is_convex_hull_edge = true;
            for(int k = 0; k < this->get_POIs_length();k++) {
                if(k != i && k != j) {
                    double c = edge.equ_substitute(this->POIs[k]);
                    if(c > edge.get_equ_c()) more++;
                    else if(c < edge.get_equ_c()) less++;
                    if(less != 0 && more != 0) {
                        is_convex_hull_edge = false;
                        break;
                    }
                }
            }
            if(is_convex_hull_edge) {
                this->add_Convex_Hull_edge(edge);
                this->POIs[i].set_isCH_POI(true);
                this->POIs[j].set_isCH_POI(true);
                break;
            }
        }
    }

    this->add_Convex_Hull_edge(Edge(this->get_CH_edge(this->get_CH_edges_length()-1).get_edge()[1], this->get_CH_edge(0).get_edge()[0]));

}

void Map::BruteForceTSP() {

    this->min_distance = INT32_MAX;
    this->shortest_route = vector<int>();
    vector<int> route;
    for(int i = 1;i < this->get_POIs_length();i++) route.push_back(i);

    do {

        vector<int> route_extend(route);
        route_extend.insert(route_extend.begin(), 0);
        route_extend.push_back(0);
        double distance = 0;
        for(int i = 0;i < route_extend.size()-1;i++) {
            distance += cal_Cartesian_distance(this->POIs[route_extend[i]], this->POIs[route_extend[i+1]]);
        }
        if(distance < this->min_distance) {
            this->min_distance = distance;
            this->shortest_route = vector<int>(route_extend);
        }

    } while(next_permutation(route.begin(), route.end()));

}

void Map::ConvexHullTSP() {

    vector<vector<POI>> append_to_which_edge(this->get_CH_edges_length());

    for(POI poi:this->POIs) {
        if(!poi.get_isCH_POI()) {
            int index = -1;
            double min_distance = INT32_MAX;
            for(int i = 0;i < this->get_CH_edges_length();i++) {
                double distance = this->Convex_Hull_edges[i].cal_distance_with_poi(poi);
                if(distance < min_distance) {
                    min_distance = distance;
                    index = i;
                }
            }
            append_to_which_edge[index].push_back(poi);
        }
    }

    for(int i = 0;i < append_to_which_edge.size();i++) sort_with_distance(append_to_which_edge[i], this->Convex_Hull_edges[i].get_edge()[0]);

    this->CH_shortest_route.push_back(this->Convex_Hull_edges[0].get_edge()[0]);
    for (int i = 0; i < append_to_which_edge.size(); i++) {
        for (int j = 0; j < append_to_which_edge[i].size(); j++) this->CH_shortest_route.push_back(append_to_which_edge[i][j]);
        this->CH_shortest_route.push_back(this->Convex_Hull_edges[i].get_edge()[1]);
    }

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
double cal_Cartesian_distance(POI poi1, POI poi2) {
    Coordinate cacoo1 = poi1.get_cartesian_coordinate(), cacoo2 = poi2.get_cartesian_coordinate();
    return sqrt(pow(cacoo1.get_longitude() - cacoo2.get_longitude(), 2) + pow(cacoo1.get_latitude() - cacoo2.get_latitude(), 2));
}
void sort_with_distance(vector<POI> &POIs, POI ref) {
    for (int i = 1;i < POIs.size();i++){
            POI key = POIs[i];
            int j = i - 1;
            while((j>=0) && (cal_Cartesian_distance(POIs[j], ref) > cal_Cartesian_distance(key, ref))) {
                    POIs[j+1] = POIs[j];
                    j--;
            }
            POIs[j+1] = key;
    }
}

