#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

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
        void set_coordinate(double, double);
        void set_coordinate(Coordinate);
        double get_longitude();
        double get_latitude();
        Coordinate get_coordinate();

    private:
        Coordinate coordinate;
};

class Map {
    public:
        Map();
        void add_POI(double, double);
        double get_POI_longitude(int);
        double get_POI_latitude(int);
        Coordinate get_POI_coordinate(int);
        int get_length();
    private:
        vector<POI> POIs;
        int length;
};

// algorithms
void BruteForceClostestPoints(Map&, Coordinate&, Coordinate&);

// tools
double cal_distance(Coordinate, Coordinate);

int main(void) {

    Map map;

    ifstream ifs;

    int number_POIs;

    ifs.open("td.txt");
    ifs >> number_POIs;

    for(int i = 0;i < number_POIs;i++) {
        double longitude, latitude;
        ifs >> longitude >> latitude;
        map.add_POI(longitude, latitude);
    }

    ifs.close();

    cout << "[POIs List]";
    for(int i = 0;i < number_POIs;i++) {
        cout << map.get_POI_longitude(i) << ' ' << map.get_POI_latitude(i) << endl;
    }
    cout << endl;

    Coordinate coo1, coo2;
    BruteForceClostestPoints(map, coo1, coo2);
    cout << "Closest Pair: (" << coo1.get_longitude() << ", " << coo1.get_latitude() << "), (" << coo2.get_longitude() << ", " << coo2.get_latitude() << ')' << endl;

    return 0;

}

// Coordinate class
Coordinate::Coordinate() { this->longitude = 0; this->latitude = 0; }
Coordinate::Coordinate(double longitude, double latitude) { set_coordinate(longitude, latitude); }
void Coordinate::set_coordinate(double longitude, double latitude) { this->longitude = longitude; this->latitude = latitude; }
double Coordinate::get_longitude() { return this->longitude; }
double Coordinate::get_latitude() { return this->latitude; }

// POI class
POI::POI() { this->coordinate = Coordinate(); }
POI::POI(double longitude, double latitude) { this->coordinate = Coordinate(longitude, latitude); }
void POI::set_coordinate(double longitude, double latitude) { this->set_coordinate(longitude, latitude); }
void POI::set_coordinate(Coordinate coordinate) { this->coordinate = coordinate; }
double POI::get_longitude() { return this->coordinate.get_longitude(); }
double POI::get_latitude() { return this->coordinate.get_latitude(); }
Coordinate POI::get_coordinate() { return this->coordinate; }

// Map class
Map::Map() { this->length = 0; }
void Map::add_POI(double longitude, double latitude) { this->POIs.push_back(POI(longitude, latitude)); this->length++; }
double Map::get_POI_longitude(int index) { return this->POIs[index].get_longitude(); }
double Map::get_POI_latitude(int index) { return this->POIs[index].get_latitude(); }
Coordinate Map::get_POI_coordinate(int index) { return this->POIs[index].get_coordinate(); }
int Map::get_length() { return this->length; }

// algorithms
void BruteForceClostestPoints(Map &map, Coordinate &coo1, Coordinate &coo2) {

    double dmin = INT32_MAX;

    for(int i = 0;i < map.get_length() - 1;i++) {
        for(int j = i+1;j < map.get_length();j++) {
            double distance = cal_distance(map.get_POI_coordinate(i), map.get_POI_coordinate(j));
            if(distance < dmin) {
                dmin = distance;
                coo1 = map.get_POI_coordinate(i);
                coo2 = map.get_POI_coordinate(j);
            }
        }
    }

}

// tools
double cal_distance(Coordinate coo1, Coordinate coo2) { return sqrt(pow(coo1.get_longitude() - coo2.get_longitude(), 2) + pow(coo1.get_latitude() - coo2.get_latitude(), 2)); }
