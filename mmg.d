import std.stdio;
import std.getopt;
import std.math;
import std.format;
import std.path;

alias M1x3 = float[3];

void m1x3_add (out M1x3 a, in M1x3 b, in M1x3 c)
{
    for (int i = 0; i < 3; i++) a[i] = b[i] + c[i];
}
void m1x3_subtract (out M1x3 a, in M1x3 b, in M1x3 c)
{
    for (int i = 0; i < 3; i++) a[i] = b[i] - c[i];
}
void m1x3_cross (out M1x3 v, in M1x3 a, in M1x3 b)
{
    v[0] = a[1]*b[2] - b[1]*a[2];
    v[1] = a[2]*b[0] - b[2]*a[0];
    v[2] = a[0]*b[1] - b[0]*a[1];
}
float m1x3_dot (in M1x3 a, in M1x3 b)
{
    float sum = a[0]*b[0];
    for (int i = 1; i < 3; i++)
    {
        sum += a[i]*b[i];
    }
    return sum;
}
bool m1x3_normalised (out M1x3 d, in M1x3 v)
{
    float m = sqrt (m1x3_dot (v, v));
    if (m <= 1e-4)
    {
        d[0] = d[1] = d[2] = 0;
        return false;
    }
    float rp = 1.0/m;
    for (int i = 0; i < 3; i++)
    {
        d[i] = rp*v[i];
    }
    return true;
}

class Point
{
    union
    {
        M1x3 xyz;
        struct
        {
            float x, y, z;
        }
    }
    this ()
    {
        x = y = z = 0;
    }
    this (float X, float Y, float Z)
    {
        x = X;
        y = Y;
        z = Z;
    }
    override size_t toHash() const @safe pure nothrow @nogc
    {
        auto X = cast (int)x&0x3ff;
        auto Y = cast (int)y&0x3ff;
        auto Z = cast (int)z&0x3ff;
        return cast (size_t)((X<<20)|(Y<<10)|(Z));
    }
    bool opEqual(ref const Point other) const @safe pure nothrow @nogc
    {
        const float epsilon = 1e-3;
        bool a = abs (x - other.x) <= epsilon;
        bool b = abs (y - other.y) <= epsilon;
        bool c = abs (z - other.z) <= epsilon;
        return a && b && c;
    }
}
class Edge
{
    public
    {
        Point vert;
        Edge next;
        Edge prev;
        Edge twin;
        Face face;
    }
    public this (Point vertex)
    {
        vert = vertex;
        next = prev = null;
        twin = null;
        face = null;
    }
}
class Face
{
    Edge loop;
    M1x3 normal;
    float distance;
    uint count;

    public this ()
    {
        loop = null;
        normal[0] = normal[1] = normal[2] = 0;
        distance = 0;
    }
    public bool append (Point vertex)
    {
        Edge e = edge_create (vertex, this);
        if (loop !is null)
        {
            loop.prev.next = e;
            e.prev = loop.prev;
            e.next = loop;
            loop.prev = e;
        }
        loop = e;
        count++;
        return true;
    }
    public void compute_plane ()
    {
        if (!valid ())
        {
            return;
        }
        Point a = loop.next.vert;
        Point b = loop.prev.vert;
        Point c = loop.next.next.vert;

        M1x3 ba, ca, x;
        m1x3_subtract (ba, b.xyz, a.xyz);
        m1x3_subtract (ca, c.xyz, a.xyz);
        m1x3_cross (x, ba, ca);
        m1x3_normalised (normal, x);
        distance = m1x3_dot (normal, a.xyz);
    }
    public bool valid () { return !(count < 3); }
}

private Point[] _points;
private Edge[] _edges;
private Edge[][Point] _edges_point;
private Face[] _faces;
private M1x3 _mins, _maxs, _centre;

Edge edge_create (Point vertex, Face face)
{
    Edge e = new Edge (vertex);
    e.prev = e.next = e;
    e.face = face;
    _edges_point[vertex] ~= e;
    _edges ~= e;
    return e;
}

void edge_connect_all ()
{
    uint connected = 0;
    uint border = 0;
    foreach (e0; _edges)
    {
        Edge next0 = e0.next;
        auto mapped = next0.vert in _edges_point;
        if (mapped is null)
        {
            continue;
        }

        foreach (e1; *mapped)
        {
            Edge next1 = e1.next;
            if (next1.vert !is e0.vert) continue;
            if (next0.vert !is e1.vert) continue;
            e1.twin = e0;
            e0.twin = e1;
            connected++;
            goto Next;
        }
        //Found a border edge
        border++;
    Next:
        continue;
    }
    writeln ("Adjacency graph:");
    writefln ("\tconnected: %s", connected);
    writefln ("\tborders: %s", border);
}

Edge[] edge_concave ()
{
    const float STEEP = 0.3;
    Edge[] result;
    foreach (e; _edges)
    {
        if (e.vert.z < _zmin)
        {
            continue;
        }
        if (e.vert.z > _zmax)
        {
            continue;
        }
        if (e.twin is null)
        {
            result ~= e;
            continue;
        }
        auto f0 = e.face;
        auto f1 = e.twin.face;
        float dot = m1x3_dot (f0.normal, f1.normal);
        if (dot >= STEEP)
        {
            continue;
        }
        result ~= e;
    }
    return result;
}

void compute_bounds ()
{
    for (int i = 0; i < 3; i++)
    {
        _mins[i] = float.infinity;
        _maxs[i] =-float.infinity;
        _centre[i] = 0;
    }
    float count = 0;
    foreach (p; _points)
    {
        for (int i = 0; i < 3; i++)
        {
            if (p.xyz[i] <= _mins[i]) _mins[i] = p.xyz[i];
            if (p.xyz[i] >= _maxs[i]) _maxs[i] = p.xyz[i];
            _centre[i] += p.xyz[i];
        }
        count += 1.0;
    }
    for (int i = 0; i < 3; i++)
    {
        _mins[i] -= _padding;
        _maxs[i] += _padding;
        _centre[i] /= count;
    }
    writefln ("Mins: %s", _mins);
    writefln ("Maxs: %s", _maxs);
    writefln ("Centroid: %s", _centre);
}

void load_bsp (string path)
{
    const int LUMP_VERTS = 2;
    const int LUMP_FACES = 6;
    const int LUMP_EDGES = 11;
    const int LUMP_SEDGE = 12;
    const int MAX_LUMPS  = 19;
    struct Lump
    {
        int offset, length;
    }
    struct Header
    {
        int id;
        int vers;
        Lump[MAX_LUMPS] lumps;
    }
    struct BSP_edge
    {
        ushort a, b;
    }
    struct BSP_face
    {
        ushort plane;
        short side;
        int edge;
        short count;
        short texinfo;
        byte[4] styles;
        int lightofs;
    }
    Header h;
    auto file = File (path, "rb");
    fread (&h, h.sizeof, 1, file.getFP ());
    
    fseek (file.getFP (), h.lumps[LUMP_VERTS].offset, SEEK_SET);
    {
        ulong count = h.lumps[LUMP_VERTS].length/M1x3.sizeof;
        for (int i = 0; i < count; i++)
        {
            M1x3 tmp;
            fread (tmp.ptr, tmp.sizeof, 1, file.getFP ());
            _points ~= new Point (tmp[0], tmp[1], tmp[2]);
        }
    }

    BSP_edge[] bedges = null;
    fseek (file.getFP (), h.lumps[LUMP_EDGES].offset, SEEK_SET);
    {
        ulong count = h.lumps[LUMP_EDGES].length/BSP_edge.sizeof;
        bedges = new BSP_edge[count];
        fread (bedges.ptr, BSP_edge.sizeof, count, file.getFP ());
    }

    int[] bsedges = null;
    fseek (file.getFP (), h.lumps[LUMP_SEDGE].offset, SEEK_SET);
    {
        ulong count = h.lumps[LUMP_SEDGE].length/int.sizeof;
        bsedges = new int[count];
        fread (bsedges.ptr, int.sizeof, count, file.getFP ());
    }

    BSP_face[] bfaces = null;
    fseek (file.getFP (), h.lumps[LUMP_FACES].offset, SEEK_SET);
    {
        ulong count = h.lumps[LUMP_FACES].length/BSP_face.sizeof;
        bfaces = new BSP_face[count];
        fread (bfaces.ptr, BSP_face.sizeof, count, file.getFP ());
    }

    foreach (ref face; bfaces)
    {
        Face f = new Face;

        //Build edge loop for the face
        int edge = face.edge;
        for (int i = 0; i < face.count; i++)
        {
            int e = bsedges[edge + i];
            if (e >= 0) f.append (_points[bedges[e].a]);
            else f.append (_points[bedges[-e].b]);
        }

        //Try to build face plane
        if (!f.valid ())
        {
            throw new Exception ("Invalid face found!");
        }
        f.compute_plane ();

        _faces ~= f;
    }

    destroy (bedges);
    destroy (bsedges);
    destroy (bfaces);
    file.close ();
}

private void write_svg (string path, Edge[] edges)
{
    auto file = File (path, "wt");
    auto width = _maxs[0] - _mins[0];
    auto height = _maxs[1] - _mins[1];
    auto hwidth = 0.5*width;
    auto hheight = 0.5*height;
    auto top = 0.0f;
    auto left = 0.0f;
    if (_raw)
    {
        top = _mins[1];
        left = _mins[0];
        width = _maxs[0];
        height = _maxs[1];
    }
    file.writefln ("<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"%s %s %s %s\">", top, left, width, height);
    foreach (e; edges)
    {
        Point a = e.vert;
        Point b = e.next.vert;

        float x0 = a.x;
        float y0 = a.y;
        float x1 = b.x;
        float y1 = b.y;
        if (!_raw)
        {
            x0 = hwidth + a.x - _centre[0];
            y0 = hheight + a.y - _centre[1];
            x1 = hwidth + b.x - _centre[0];
            y1 = hheight + b.y - _centre[1];
        }

        file.writefln (
            "<line x1=\"%s\" y1=\"%s\" x2=\"%s\" y2=\"%s\" stroke=\"black\"/>",
            x0, y0,
            x1, y1
        );
    }
    file.writefln ("</svg>");
    file.close ();
}

private float _padding = 256.0;
private float _zmax = float.infinity;
private float _zmin =-float.infinity;
private bool _raw = false;
private string _prefix = "";

void main (string[] args)
{
    auto info = getopt (
        args,
        std.getopt.config.passThrough,
        "padding", "Padding around the origin for an edge to be printed", &_padding,
        "zmin", "Lower bound of vertices for consideration", &_zmin,
        "zmax", "Upper bound of vertices for consideration", &_zmax,
        "raw", "Keep lines centered around origin", &_raw,
        "prefix", "Prefix to append to output", &_prefix
    );
    if (info.helpWanted || args.length < 2)
    {
        defaultGetoptPrinter ("USAGE: mmg <opts> map", info.options);
        return;
    }

    load_bsp (args[$-1]);
    compute_bounds ();
    
    edge_connect_all ();
    
    auto edges = edge_concave ();
    writefln ("Concave edges: %s", edges.length);

    write_svg (_prefix ~ baseName (stripExtension (args[$-1])) ~ ".svg", edges);
}
