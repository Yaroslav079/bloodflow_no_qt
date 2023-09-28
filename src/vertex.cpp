#include "vertex.h"

Meta_vertex::Meta_vertex(const std::string & id)
    : id(id)
{}
Meta_vertex::~Meta_vertex()
{}

std::string Meta_vertex::get_id()
{
    return id;
}
