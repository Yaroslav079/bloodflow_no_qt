#include "graph.h"
#include <cassert>
#include <iostream>

Simple_vertex::Simple_vertex(const std::string & id, double x, double y, double z)
    : id(id), x(x), y(y), z(z)
{}
double Simple_vertex::get_x() const
{
    return x;
}

double Simple_vertex::get_y() const
{
    return y;
}

double Simple_vertex::get_z() const
{
    return z;
}

void Simple_vertex::set_x(double value)
{
    x = value;
}

void Simple_vertex::set_y(double value)
{
    y = value;
}

void Simple_vertex::set_z(double value)
{
    z = value;
}

std::string Simple_vertex::get_id() const
{
    return id;
}

std::ostream & operator<<(std::ostream &os, const Simple_vertex &sv)
{
    os << "Vertex " << sv.get_id();
    return os;
}


Simple_vertex *Simple_edge::get_first_vertex()
{
    return first_vertex;
}
Simple_vertex *Simple_edge::get_last_vertex()
{
    return last_vertex;
}
Simple_edge::Simple_edge(const std::string & id, Simple_vertex *first,
                                    Simple_vertex *last)
:id(id), first_vertex(first), last_vertex(last)
{
    assert(first_vertex != nullptr && last_vertex != nullptr && first_vertex != last_vertex);
}

void Simple_edge::set_vertices(Simple_vertex *first, Simple_vertex *last)
{
    first_vertex = first;
    last_vertex = last;
}

std::ostream & operator<<(std::ostream &os, const Simple_edge &se)
{
    os << "Edge " << se.get_id() << " made of " << *se.first_vertex << " and "
              << *se.last_vertex;
    return os;
}

std::string Simple_edge::get_id() const
{
    return id;
}
