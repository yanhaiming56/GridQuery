#ifndef INECLASSTYPE_H
#define INECLASSTYPE_H

#include "INEBeautyCode.h"

class INEPoint {
public:
    INEPoint() : m_X(-1), m_Y(-1) {}

    INEPoint(double x, double y) : m_X(x), m_Y(y) {}

    INEPoint(const INEPoint &point) {
        m_X = point.m_X;
        m_Y = point.m_Y;
    }

    INEPoint &operator=(INEPoint &point) {
        m_X = point.m_X;
        m_Y = point.m_Y;
        return *this;
    }

    double m_X;
    double m_Y;
};

class INEVertex : public INEPoint {
public:
    INEVertex() : m_Vid(-1) {}

    INEVertex(int id, double x, double y) : INEPoint(x, y) {
        m_Vid = id;
    }

    INEVertex(const INEVertex &v) {
        m_X = v.m_X;
        m_Y = v.m_Y;
        m_Vid = v.m_Vid;
        m_Roads = v.m_Roads;
    }

    INEVertex &operator=(const INEVertex &v) {
        m_X = v.m_X;
        m_Y = v.m_Y;
        m_Vid = v.m_Vid;
        m_Roads = v.m_Roads;
        return *this;
    }

    int m_Vid; //顶点编号
    set<int> m_Roads; //The roads cover the vertex

};

class Entity : public INEPoint {
public:
    Entity() : m_Oid(-1) {}

    Entity(int id, double x, double y) : INEPoint(x, y) {}

    Entity(const Entity &entity) {
        m_X = entity.m_X;
        m_Y = entity.m_Y;
        m_Oid = entity.m_Oid;
    }

    Entity &operator=(const Entity &entity) {
        m_X = entity.m_X;
        m_Y = entity.m_Y;
        m_Oid = entity.m_Oid;
        return *this;
    }

    int m_Oid; //实体点编号
};

class INERoad {
public:
    INERoad() : m_Rid(-1) {}

    INERoad(int id, int startId, int endId, double distance) : m_Rid(id), m_StartId(startId), m_EndId(endId),
                                                            m_Distance(distance) {}

    INERoad(const INERoad &road) {
        m_Rid = road.m_Rid;
        m_StartId = road.m_StartId;
        m_EndId = road.m_EndId;
        m_Distance = road.m_Distance;
        m_Entitys = road.m_Entitys;
    }

    INERoad &operator=(const INERoad &road) {
        m_Rid = road.m_Rid;
        m_StartId = road.m_StartId;
        m_EndId = road.m_EndId;
        m_Distance = road.m_Distance;
        m_Entitys = road.m_Entitys;
        return *this;
    }

    int m_Rid; //道路编号
    int m_StartId; //起点
    int m_EndId; //终点
    double m_Distance; //道路长度
    set<int> m_Entitys; //道路上实体点

};


class Data {
public:
    Data() : m_Oid(-1), m_Distance(-1) {}

    Data(int id, double distance) : m_Oid(id), m_Distance(distance) {}

    Data(const Data &data) {
        m_Oid = data.m_Oid;
        m_Distance = data.m_Distance;
    }

    Data &operator=(const Data &data) {
        m_Oid = data.m_Oid;
        m_Distance = data.m_Distance;
        return *this;
    }

    friend bool operator<(const Data &data1, const Data &data2) {
        if (data1.m_Oid == data2.m_Oid)
            return false;

        if (data1.m_Distance < data2.m_Distance)
            return true;

        else
            return false;
    }

    int m_Oid; //实体点ID
    double m_Distance; //距离起始顶点距离
};

class VisitedNode {
public:
    VisitedNode() : m_Vid(-1), m_Distance(-1) {}

    VisitedNode(int id, double distance) : m_Vid(id), m_Distance(distance) {}

    VisitedNode(const VisitedNode &vn) {
        m_Vid = vn.m_Vid;
        m_Distance = vn.m_Distance;
    }

    VisitedNode &operator=(const VisitedNode &vn) {
        m_Vid = vn.m_Vid;
        m_Distance = vn.m_Distance;
        return *this;
    }

    friend bool operator<(const VisitedNode &vn1, const VisitedNode &vn2) {
        if (vn1.m_Vid == vn2.m_Vid)
            return false;

        if (vn1.m_Distance > vn2.m_Distance)
            return true;

        else
            return false;
    }

    int m_Vid; //已访问的节点
    double m_Distance; //查询点到节点距离
};


class VisitedEntity {
public:
    VisitedEntity() : m_Oid(-1), m_Distance(-1) {}

    VisitedEntity(int id, double distance) : m_Oid(id), m_Distance(distance) {}

    VisitedEntity(const VisitedEntity &ve) {
        m_Oid = ve.m_Oid;
        m_Distance = ve.m_Distance;
    }

    VisitedEntity &operator=(const VisitedEntity &ve) {
        m_Oid = ve.m_Oid;
        m_Distance = ve.m_Distance;
        return *this;
    }

    friend bool operator<(const VisitedEntity &ve1, const VisitedEntity &ve2) {
        if (ve1.m_Oid == ve2.m_Oid)
            return false;
        if (ve1.m_Distance < ve2.m_Distance)
            return true;
        else
            return false;
    }

    int m_Oid;
    double m_Distance; //The distance from query point
};

#endif // CLASSTYPE_H
