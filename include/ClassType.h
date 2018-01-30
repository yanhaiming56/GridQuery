#ifndef CLASSTYPE_H
#define CLASSTYPE_H

#include "CodeBeauty.h"


using namespace std;

struct Point
{
    double m_XPos;
    double m_YPos;

    Point()
    {
        m_XPos = -1;
        m_YPos = -1;
    }

    Point( double xPos, double yPos )
    {
        m_XPos = xPos;
        m_YPos = yPos;
    }

    Point( const Point& point )
    {
        m_XPos = point.m_XPos;
        m_YPos = point.m_YPos;
    }

    friend ostream& operator << ( ostream& out, const Point& point )
    {
        out << point.m_XPos << " " << point.m_YPos;
        return out;
    }
};

struct Vertex
{
    int m_Vid;
    int m_Oid;
    vector<int> m_Rids;
    Point m_Point;

    Vertex()
    {
        m_Vid = -1;
        m_Oid = -1;
        m_Point = Point();
    }

    Vertex( const Vertex& vertex )
    {
        m_Vid = vertex.m_Vid;
        m_Oid = vertex.m_Oid;
        m_Rids = vertex.m_Rids;
        m_Point = vertex.m_Point;
    }

    Vertex( int vid, Point point ) : m_Vid( vid ), m_Oid( -1 ), m_Point( point )
    {
    }

    friend ostream& operator << ( ostream& out, const Vertex& vertex )
    {
        out << vertex.m_Vid << " " << vertex.m_Point;
        return out;
    }
};

struct DataPoint
{
    int m_Oid;
    double m_VDistance;
    Point m_Point;

    DataPoint()
    {
        m_Oid = -1;
        m_VDistance = -1;
        m_Point = Point();
    }
    DataPoint( int oid, double distance, Point point ) : m_Oid( oid ), m_VDistance( distance ),
        m_Point( point ) {}

    DataPoint( const DataPoint& dataPoint )
    {
        m_Oid = dataPoint.m_Oid;
        m_VDistance = dataPoint.m_VDistance;
        m_Point = dataPoint.m_Point;
    }

    bool operator < ( const DataPoint& dataPoint ) const
    {
        if ( m_VDistance < dataPoint.m_VDistance )
            return true;

        else
            return false;
    }
    friend ostream& operator << ( ostream& out, const DataPoint& dataPoint )
    {
        out << dataPoint.m_Oid << " " << dataPoint.m_Point;
        return out;
    }
};

struct Road
{
    int m_Rid;
    double m_Distance;
    double m_Slope;
    Vertex m_Begin;
    Vertex m_End;
    set<DataPoint> m_Datas;

    Road()
    {
        m_Rid = -1;
        m_Distance = -1;
        m_Begin = Vertex();
        m_End = Vertex();
    }

    Road( const Road& road )
    {
        m_Rid = road.m_Rid;
        m_Begin = road.m_Begin;
        m_End = road.m_End;
        m_Datas = road.m_Datas;
        m_Slope = road.m_Slope;
        m_Distance = road.m_Distance;
    }
    Road( int rid, double distance, Vertex startPoint, Vertex endPoint )
    {
        m_Rid = rid;
        m_Distance = distance;
        m_Begin = startPoint;
        m_End = endPoint;
        if ( fabs( startPoint.m_Point.m_XPos - endPoint.m_Point.m_XPos ) < ZERO )
            m_Slope = MAXNUM;
        else
            m_Slope = ( startPoint.m_Point.m_YPos - endPoint.m_Point.m_YPos ) / ( startPoint.m_Point.m_XPos -
                      endPoint.m_Point.m_XPos );
    }

    friend ostream& operator << ( ostream& out, const Road& road )
    {
        out << road.m_Rid << " " << road.m_Begin.m_Vid << " " << road.m_End.m_Vid << " " << road.m_Distance << " " << road.m_Datas.size();
        return out;
    }
};

struct SubRoad
{
    int m_Rid;
    int m_Oid;
    Point m_Start;
    Point m_End;
    double m_Slope;
    bool m_IsVertex;

    SubRoad() : m_Rid( -1 ), m_Oid( -1 ), m_Slope( 0 ), m_IsVertex( false ) {}

    SubRoad( int rid, int oid, Point startPoint, Point endPoint, double slope,
             bool isVertex ) : m_Rid( rid ),
        m_Oid( oid ), m_Slope( slope ), m_IsVertex( isVertex )
    {
        m_Start = startPoint;
        m_End = endPoint;
    }

    SubRoad( const SubRoad& subRoad )
    {
        m_Rid = subRoad.m_Rid;
        m_Oid = subRoad.m_Oid;
        m_Slope = subRoad.m_Slope;
        m_Start = subRoad.m_Start;
        m_End = subRoad.m_End;
        m_IsVertex = subRoad.m_IsVertex;

    }
};

struct Grid
{
    int m_X;
    int m_Y;
    vector<SubRoad> subRoads;

    Grid() : m_X( -1 ), m_Y( -1 ) {}
    Grid( int x, int y ) : m_X( x ), m_Y( y ) {}
    Grid( const Grid& grid )
    {
        m_X = grid.m_X;
        m_Y = grid.m_Y;
        subRoads.insert( subRoads.end(), grid.subRoads.begin(), grid.subRoads.end() );
    }
    Grid operator = ( const Grid& grid )
    {
        m_X = grid.m_X;
        m_Y = grid.m_Y;
        subRoads.insert( subRoads.end(), grid.subRoads.begin(), grid.subRoads.end() );
        return *this;
    }
    void insertSubRoad( const int& x, const int& y, const SubRoad& subRoad )
    {
        m_X = x;
        m_Y = y;
        subRoads.push_back( subRoad );
    }
};
struct Vertexs
{
    Vertex* m_Vertexs;
    Vertexs()
    {
        m_Vertexs = nullptr;
    }
};

struct DataPoints
{
    DataPoint* m_DataPoints;
    DataPoints()
    {
        m_DataPoints = nullptr;
    }
};

struct Roads
{
    Road* m_Roads;
    Roads()
    {
        m_Roads = nullptr;
    }
};

struct SubRoads
{
    SubRoad* m_SubRoads;
    SubRoads()
    {
        m_SubRoads = nullptr;
    }
};

struct Grids
{
    Grid* m_Grids;
    Grids()
    {
        m_Grids = nullptr;
    }
};


#endif // CLASSTYPE_H
