#include "edge.h"

Meta_edge::Meta_edge(const std::string & name, Simple_vertex *first, Simple_vertex *last)
    : Simple_edge(name, first, last)
{}
Meta_edge::~Meta_edge()
{}
