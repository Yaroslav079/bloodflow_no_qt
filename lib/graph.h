#ifndef GRAPH_H
#define GRAPH_H
#include <string>


#define DUMP_TO_OS(x) os << (#x) << ": " << (x) << std::endl

class Simple_edge;

class Simple_vertex
{
    /*!
     * \brief A base class representing a vertex of a CV graph only as a geometric entity.
     */
    std::string id;
    double x, y, z;



public:
    Simple_vertex(const std::string & id, double x = 0, double y = 0, double z = 0);
    double get_x() const;
    double get_y() const;
    double get_z() const;
    void set_x(double value);
    void set_y(double value);
    void set_z(double value);
    std::string get_id() const;
    friend std::ostream& operator<<(std::ostream& os, const Simple_vertex& sv);
};

typedef enum {flow, speed, area, pressure} quantity; /*!< to obtain different physical quantities*/

class Simple_edge
{
    /*!
     * \brief A base class representing an edge of a CV graph.
     */
    std::string id;
protected:
    Simple_vertex *first_vertex;
    Simple_vertex *last_vertex;
public:
    Simple_vertex* get_first_vertex();
    Simple_vertex* get_last_vertex();
    Simple_edge(const std::string & id, Simple_vertex *first = nullptr,
                                        Simple_vertex *last = nullptr);
    void set_vertices(Simple_vertex *first = nullptr, Simple_vertex *last = nullptr);
    std::string get_id() const;
    friend std::ostream& operator<<(std::ostream& os, const Simple_edge& se);

    virtual const double * get_data(quantity vt, int & size, int & step) = 0;
    virtual double get_P_unsafe(int i) = 0;
    virtual double get_Q_unsafe(int i) const = 0;
    virtual double get_S_unsafe(int i) const = 0;
    virtual double get_U_unsafe(int i) const = 0;
    virtual int get_points_num() const = 0;

};
#endif
