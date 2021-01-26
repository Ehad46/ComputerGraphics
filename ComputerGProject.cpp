#include<iostream>
#include<vector>
#include<math.h>
#include<fstream>
#include<sstream>
using namespace std;

class Material
{
    // valores = { ambiente, difuso, especular, exponente }
    vector<double> valores;
    public: 
    Material(double a, double d,double s,double ns)
    {
        valores.push_back(a);
        valores.push_back(d);
        valores.push_back(s);
        valores.push_back(ns);
    }
    vector<double> obtenerValores()
    {
        return valores;
    } 
};
class Luz
{
    vector<double> position;
    vector<int> color;
    vector<double> attenVar; // 0, 1 y 2
    vector<double> direction; // En caso de usar cono
    double angle; // Angulo de apertura
    bool isCone = false; 
    public: 
    Luz(vector<double> p, vector<int> c, vector<double> v)
    {
        position = p; color = c; attenVar = v; isCone = false;
    }
    Luz(vector<double> p, vector<int> c, vector<double> v, vector<double> dir, double a)
    {
        position = p; color = c; attenVar = v; direction = dir; isCone = true; angle = a;
    }
    vector<double> regresarPosicion() { return position; }
    vector<int> regresarColor() { return color; }
    vector<double> regresarAtenuacion() { return attenVar; }
    vector<double> regresarDireccion() { return direction; }
    double regresarAngulo() { return angle; }
    bool is_Cone() { return isCone; }
};
class Transformacion
{
    vector<vector<double> > tabla;
    public:
    Transformacion(float f){ // Transformacion de escala en todos los ejes
        tabla = vector<vector<double> >(0);
        for (int i = 0; i < 4; i++)
        {
            vector<double> subtabla;
            for (int i2 = 0; i2 < 4; i2++)
                subtabla.push_back((i==i2)? ((i!=3)? f : 1) : 0); // Si estoy en la diagonal identidad, meter f, sino 0
            tabla.push_back(subtabla);
        }           
        /*
        fx 0 0 0
        0 fy 0 0
        0 0 fz 0
        0 0 0 1
        */
    }
    Transformacion(float x, float y, float z){ // Tres numeros , entonces es una transformacion de movimiento
        tabla = vector<vector<double> >(4, vector<double>(4, 0));
        for (int i = 0; i < 4; i++)
            tabla[i][i] = 1;
        tabla[0][3] = x;
        tabla[1][3] = y;
        tabla[2][3] = z;
    }
    Transformacion(char axis, float angle){ // Un eje y un numero, entonces es una transformacion de rotacion
        tabla = vector<vector<double> >(4, vector<double>(4, 0));
        for (int i = 0; i < 4; i++)
            tabla[i][i] = 1;
        switch (axis)
        {
            case 'x':
                tabla[1][1] = cos(angle);
                tabla[2][1] = sin(angle);
                tabla[1][2] = - tabla[2][1];
                tabla[2][2] = tabla[1][1];
                break;
            case 'y':
                tabla[0][0] = cos(angle);
                tabla[0][2] = sin(angle);
                tabla[2][0] = - tabla[0][2];
                tabla[2][2] = tabla[0][0];
                break;
            default: // Entonces es z 
                tabla[0][0] = cos(angle);
                tabla[1][0] = sin(angle);
                tabla[0][1] = - tabla[1][0];
                tabla[1][1] = tabla[0][0]; 
                break;
        }     
    }
    vector<vector<double> > regresar_Tabla() { return tabla; }
    void definir_Tabla(vector<vector<double> > nueva) { tabla = nueva; }
};
class OBJ
{
    string nombreObjeto = " ";
    vector<vector<double> > puntos;
    vector<vector<double> > caras;
    vector<vector<double> > lados;
    vector<vector<double> > normales;
    vector<vector<double> > normalesVertices;
    Material mt = Material(1,0.5,0.5,150);
    vector<int> rgb = vector<int>(3,0);

    public:
    OBJ(string string, bool i){ // Nombre y bool para determinar si es vlf u obj (true - vlf)
      nombreObjeto = string;
      if(i) 
      {
          lecturaVlf(string); 
      }
      else 
      {
          lecturaObj(string); 
      }
    }

    // Para obtener las cosas ya que son privadas
    vector<vector<double> > regresarPuntos(){ 
        return puntos; 
    }
    vector<vector<double> > regresarLados() { 
        return lados; 
    }
    vector<vector<double> > regresarCaras() { 
        return caras; 
    }
    vector<vector<double> > regresarNormales() { 
        return normales; 
    }
    vector<int> regresarRGB() { 
        return rgb; 
    }
    Material regresarMaterial() { 
        return mt; 
    }
    vector<vector<double> > regresarNormalesVertices() { 
        return normalesVertices; 
    }

    void cambiar_color(vector<int> c){ 
        for(int i=0;i<3;i++)
        {
            if(c[i]>255) rgb[i] = 255;
            else if(c[i]<0)   rgb[i] = 0;
            else rgb[i] = c[i];
        }
    }    
    void cambiar_Material(Material material) { 
        mt = material; 
    }
    bool lecturaObj(string s){
        ifstream archivo(s + ".obj");
        if (archivo.is_open()){
            string linea;
            while (getline(archivo, linea)){
                istringstream indv(linea);
                string num;
                getline(indv, num, ' ');
                string type = num;
                if (type == "v"){
                    vector<double> xyzw;
                    while (getline(indv, num, ' '))
                        xyzw.push_back(stof(num));
                    puntos.push_back(xyzw);
                }
                if (type == "f"){
                    vector<double> faces;
                    while (getline(indv, num, ' '))
                        faces.push_back((double) stoi(num));
                    caras.push_back(faces);
                }
            }
            for(vector<double> &p : puntos){
                if(p.size() < 4) p.push_back(1);
                p[1] *= -1;
            }
        }
        else
        { cout << "Error!"; return false; }
        archivo.close();
        return true;
    }
    bool lecturaVlf(string s){
        ifstream archivo(s + ".vlf", ios::binary | ios::in);
        if (archivo.is_open()){ // # puntos, # lados, # caras
            int tamPuntos, tamLados, tamCaras;
            archivo.read((char *)&tamPuntos, sizeof(int));
            archivo.read((char *)&tamLados, sizeof(int));
            archivo.read((char *)&tamCaras, sizeof(int));
            vector<vector<double> > new_puntos(tamPuntos, vector<double>(3, 0));
            vector<vector<double> > new_lados(tamLados, vector<double>(2, 0));
            vector<vector<double> > new_caras(tamCaras, vector<double>(3, 0));
            vector<vector<double> > new_normales(tamCaras, vector<double>(3, 0));
            vector<vector<double> > new_vertexNormals(tamPuntos, vector<double>(3, 0));
            for (int i = 0; i < tamPuntos; i++)
                for (int cnt = 0; cnt < 3; cnt++)
                    archivo.read((char *)&new_puntos[i][cnt], sizeof(double));
            for (int i = 0; i < tamLados; i++)
                for (int cnt = 0; cnt < 2; cnt++)
                    archivo.read((char *)&new_lados[i][cnt], sizeof(double));
            for (int i = 0; i < tamCaras; i++)
                for (int cnt = 0; cnt < 3; cnt++)
                    archivo.read((char *)&new_caras[i][cnt], sizeof(double));
            for (int i = 0; i < tamCaras; i++)
                for (int cnt = 0; cnt < 3; cnt++)
                    archivo.read((char *)&new_normales[i][cnt], sizeof(double));
            for (int i = 0; i < tamPuntos; i++)
                for (int cnt = 0; cnt < 3; cnt++)
                    archivo.read((char *)&new_vertexNormals[i][cnt], sizeof(double));
            puntos = new_puntos;
            lados = new_lados;
            caras = new_caras;
            normales = new_normales;
            normalesVertices = new_vertexNormals;
            for (int i = 0; i < puntos.size(); i++) 
                puntos[i].push_back(1);
            for (int i = 0; i < normales.size(); i++)
                normales[i].push_back(1);
            for (int i = 0; i < normalesVertices.size(); i++)
                normalesVertices[i].push_back(1);
        }
        else
        {
            cout << "Error!";
            return false;
        }
        archivo.close();
        return true;
    }
    void aplicar_Transformacion(Transformacion trns)
    {
        vector<vector<double> > t = trns.regresar_Tabla();
        for (vector<double> &p : puntos)
        {
            vector<double> nuevop(4,0);
            for (int i = 0; i < 4; i++)
            {
                float valor = 0;
                for (int i2 = 0; i2 < 4; i2++)
                {
                    valor += t[i][i2] * p[i2];
                }
                nuevop[i] = valor;
            }
            p = nuevop;
        }
    }
    void centrar(int rx, int ry){ // Para camara ortogonal
        float xmin = 0;
        float ymin = 0;
        float xmax = 0;
        float ymax = 0;
        for (vector<double> p : puntos){ // Para determinar cuanto se escala
            if (p[0] < xmin)
                xmin = p[0];
            if (p[1] < ymin)
                ymin = p[1];

            if (p[0] > xmax)
                xmax = p[0];
            if (p[1] > ymax)
                ymax = p[1];
        }
        float dx = xmax - xmin;
        float dy = ymax - ymin;
        float factor = (rx<ry)? rx/dx : ry/dy;
        Transformacion escalar(factor);
        aplicar_Transformacion(escalar); // Aplicar escala
        ymin = 0;
        xmin = 0;
        float zmin = 0;
        for (vector<double> p : puntos){ // Se busca la x y la y minima
            if (p[0] < xmin)
                xmin = p[0];
            if (p[1] < ymin)
                ymin = p[1];
            if (p[2] < zmin)
                zmin = p[2];
        }
        Transformacion mover(-xmin, -ymin, -zmin); // Se pregunta si alguno delos ejes es negativo, en caso de serlo se mueve
        aplicar_Transformacion(mover);
    }
    void dividir_Homogeneo(){ // Salir del homogeneous space
        for (vector<double> &p : puntos)
            for (int i = 0; i < 3; i++)
                p[i] /= p[3];
    }
    void triangular(){
        for (int ci = 0; ci < caras.size(); ci++)
        {
            vector<double> caraAnalizada = caras[ci];
            if (caraAnalizada.size() > 3) // Si tiene más de 3 lados
            {
                int cnt = 2;
                for (int i = caraAnalizada.size() - 3; i > 0; i--) // Se determina el numero de caras nuevas y se toman los puntos de la original
                {
                    vector<double> nueva_cara(3, 0);
                    nueva_cara[0] = caraAnalizada[0];
                    nueva_cara[1] = caraAnalizada[cnt];
                    nueva_cara[2] = caraAnalizada[cnt + 1];
                    cnt++;
                    caras.push_back(nueva_cara);
                }
            }
            caras[ci].resize(3);
        }
    }
    void calcular_Normales(){
        normales.clear();
        for(vector<double> c : caras)
        {
            vector<double> normal(4,1);        
            double v1x, v1y, v1z, v2x, v2y, v2z;
            // A cada edge se le saca sus dos puntos y le resta el primero al segundo
            v1x = puntos[ lados[ c[0] - 1 ][1] - 1 ][0] - puntos[ lados[ c[0] - 1 ][0] - 1 ][0];
            v1y = puntos[ lados[ c[0] - 1 ][1] - 1 ][1] - puntos[ lados[ c[0] - 1 ][0] - 1 ][1];
            v1z = puntos[ lados[ c[0] - 1 ][1] - 1 ][2] - puntos[ lados[ c[0] - 1 ][0] - 1 ][2];
            v2x = puntos[ lados[ c[1] - 1 ][1] - 1 ][0] - puntos[ lados[ c[1] - 1 ][0] - 1 ][0];
            v2y = puntos[ lados[ c[1] - 1 ][1] - 1 ][1] - puntos[ lados[ c[1] - 1 ][0] - 1 ][1];
            v2z = puntos[ lados[ c[1] - 1 ][1] - 1 ][2] - puntos[ lados[ c[1] - 1 ][0] - 1 ][2];
            // Calculo de vectores
            normal[0] = -1 * ((v1y * v2z) - (v1z * v2y));
            normal[1] = -1 * ((v1z * v2x) - (v1x * v2z));
            normal[2] = -1 * ((v1x * v2y) - (v1y * v2x));

            normales.push_back(normal);
        }
    }
    void calcular_NormalesVertices(){
        normalesVertices.clear(); // Se limpia para recalcular
        vector<vector<double> > normalesOrigen = regresarNormales();
        for(int i=0; i < puntos.size(); i++)
        {
            vector<double> vn(3,0);
            vector<vector<double> > vadj; // vectores adyacentes
            for(int c=0; c < caras.size(); c++) // Se busca lados a los que está conectado el punto
            {
                bool includes = false;
                for(int l=0; l<2; l++)
                {
                    if( (lados[caras[c][l] - 1][0] - 1) ==  i ) includes = true;
                    if( (lados[caras[c][l] - 1][1] - 1) ==  i ) includes = true;
                }
                if(includes) vadj.push_back(normalesOrigen[c]);
            }   
            for(vector<double> v : vadj) // La normal del punto es un promedio del valor de los vectores de los vectores adyacentes
            {
                for (int i = 0; i < 3; i++)
                    vn[i] = vn[i] + v[i];
            }
            for(int cnt=0; cnt<3; cnt++)
                vn[cnt] /= vadj.size();
            normalesVertices.push_back(vn);
        }
    }
    void transformar_VLF(){ // Incluye la triangulacion , las normales, y las normales en los vertices
        triangular();
        vector<vector<double> > nuevasCaras;
        for (int ci = 0; ci < caras.size(); ci++) // Cambio de formato a VLF, las caras usaran lados en lugar de puntos
        {
            int carastam = caras.size();
            vector<double> face = caras[ci];
            vector<double> newFace;
            for (int cnt = 0; cnt < 3; cnt++)
            {
                vector<double> edge(2, 0);
                edge[0] = face[cnt];
                edge[1] = face[(cnt + 1) % 3];
                lados.push_back(edge);
                newFace.push_back(lados.size());
            }
            nuevasCaras.push_back(newFace);
        }
        caras = nuevasCaras;
        calcular_Normales();
        calcular_NormalesVertices();
        ofstream archivo(nombreObjeto + ".vlf", ios::binary | ios::out);
        // Crear un archivo con el mismo nombre del objeto
        if (archivo.is_open())
        {
            int tam = puntos.size();
            archivo.write((char *)&tam, sizeof(int));
            tam = lados.size();
            archivo.write((char *)&tam, sizeof(int));
            tam = caras.size();
            archivo.write((char *)&tam, sizeof(int));

            for (int i = 0; i < puntos.size(); i++)
                for (int cnt = 0; cnt < 3; cnt++)
                    archivo.write((char *)&puntos[i][cnt], sizeof(double));
            for (int i = 0; i < lados.size(); i++)
                for (int cnt = 0; cnt < 2; cnt++)
                    archivo.write((char *)&lados[i][cnt], sizeof(double));
            for (int i = 0; i < caras.size(); i++)
                for (int cnt = 0; cnt < 3; cnt++)
                    archivo.write((char *)&caras[i][cnt], sizeof(double));
            for (int i = 0; i < normales.size(); i++)
                for (int cnt = 0; cnt < 3; cnt++)
                    archivo.write((char *)&normales[i][cnt], sizeof(double));
            for (int i = 0; i < puntos.size(); i++)
                for (int cnt = 0; cnt < 3; cnt++)
                    archivo.write((char *)&normalesVertices[i][cnt], sizeof(double));
        }
        else
        {
            cout << "Error!";
        }
        archivo.close();
    }
    void recalcular_Normales(){
        calcular_Normales();
        calcular_NormalesVertices();
    }
};
class operacionMatematica
{
    public:
    operacionMatematica(){ }
    double magnitud_Vector(vector<double> u)
    {
        return (sqrt(pow(u[0], 2) + pow(u[1], 2) + pow(u[2], 2)));
    }
    double productoPunto(vector<double> u, vector<double> v)
    {
        double ans = 0;
        for (int i = 0; i < 3; i++)
        {
            ans += u[i] * v[i];
        }
        return ans;
    }
    vector<double> restarVector(vector<double> u, vector<double> v)
    {
        vector<double> nuevo(4, 1);
        for (int i = 0; i < 3; i++)
            nuevo[i] = u[i] - v[i];
        return nuevo;
    }
    vector<double> reflejoVector(vector<double> u, vector<double> n)
    {
        vector<double> r(3,0);
        double k = 2 * productoPunto(u,n) * (1 / (sqrt(pow(n[0], 2) + pow(n[1], 2) + pow(n[2], 2))));
        for(int i = 0; i < 3; i++)
            r[i] = u[i] - (k * n[i]);
        return r;
    }
};
operacionMatematica op = operacionMatematica();
class Raster
{
    vector<vector<vector<unsigned char> > > image;
    vector<vector<vector<double> > > xyPixel; // Guarda la posicion x, y del modelo original
    vector<vector<double> > zBuffer;
    vector<vector<Material> > materialPixel;
    vector<vector<vector<double> > > normalPixel;
    int rx, ry;
    unsigned char r, g , b;
    float f = 500; // La camara es en perspectiva y este es el focal lenght
    vector<int> luzAmbiental = vector<int>(3,0);
    
    public:
    Raster()
    {

    }
    Raster(int _rx, int _ry){
        rx=_rx; ry=_ry;
        image = vector<vector<vector<unsigned char> > >(rx, vector<vector<unsigned char> >(ry, vector<unsigned char>(3, 0x00)));
        xyPixel = vector<vector<vector<double> > >(rx, vector<vector<double> >(ry, vector<double>(2, (double) INFINITY)));
        normalPixel = vector<vector<vector<double> > >(rx, vector<vector<double> >(ry, vector<double>(3, 0)));
        materialPixel = vector<vector<Material> >(rx,vector<Material>(ry, Material(0,0,0,0)) );
        zBuffer = vector<vector<double> >(rx, vector<double>(ry,  (double) INFINITY ));
    }
    void setTamano_Imagen(int _rx, int _ry)
    {
        rx = _rx;
        ry = _ry;
    }
    void crearPPM(string s)
    {
        ofstream imagen(s + ".ppm");
        imagen << "P6\n"
            << to_string(rx) << " " << to_string(ry) << " 255\n";
        for (int iy = 0; iy < ry; iy++)
            for (int ix = 0; ix < rx; ix++)
                imagen << image[ix][iy][0] << image[ix][iy][1] << image[ix][iy][2];
        imagen.close();
    }

    void setTamano_Focal(float foc){ 
        f = foc; 
    }
    void setColor(int red, int green, int blue){   
        r = red;     
        g = green;     
        b = blue;  
    }
    void pixel(int x, int y)
    {  
        if(x>=rx||y>=ry||y<0||x<0) 
            return;
        image[x][y][0] = r;
        image[x][y][1] = g;
        image[x][y][2] = b;
    }
    // Practica 6 : ZBuffer
    void pixel(int x, int y, double z, vector<int> color) 
    {  
        if(x>=rx||y>=ry||y<0||x<0) 
            return;
        if(z < zBuffer[x][y])
        {
            image[x][y][0] = color[0];
            image[x][y][1] = color[1];
            image[x][y][2] = color[2];
            zBuffer[x][y] = z;
        }
    }
    // Guarda el valor z y tambien los valores en x y originales del OBJ
    void pixel(int ix, int iy, int ox, int oy, double z)
    {  
        if(ix>=rx||iy>=ry||iy<0||ix<0) return;
        if(z < zBuffer[ix][iy])
        {
            image[ix][iy][0] = r;
            image[ix][iy][1] = g;
            image[ix][iy][2] = b;
            zBuffer[ix][iy] = z;
            // Actualiza valor x, y original
            xyPixel[ix][iy][0] = ox;
            xyPixel[ix][iy][1] = oy;
        }
    }
    // Guarda z, x y original y ademas tambien lleva registro del material y normal
    void pixel(int ix, int iy, int ox, int oy, double z, vector<int> color, Material mt, vector<double> n)
    {
        if (ix >= rx || iy >= ry || iy < 0 || ix < 0)
            return;
        if (z < zBuffer[ix][iy])
        {
            image[ix][iy][0] = color[0];
            image[ix][iy][1] = color[1];
            image[ix][iy][2] = color[2];
            zBuffer[ix][iy] = z;
            // Actualiza valor x, y original
            xyPixel[ix][iy][0] = ox;
            xyPixel[ix][iy][1] = oy;
            // Material y normal del pixel
            materialPixel[ix][iy] = mt;
            normalPixel[ix][iy] = n;
        }
    }
    
    // Lineas
    void linea_Naive(int x1, int y1, int x2, int y2) // Practica 1
    {
        pixel(x1, y1);
        if (x1 == x2 && y1 == x2) return;
        pixel(x2, y2);
        float m, dx, dy;
        dy = y2 - y1;
        dx = x2 - x1;
        if (abs((int)dx) >= abs((int)dy))// Si X > Y, Y por cada X
        {
            if (x1 > x2)
            {
                swap(y1, y2);
                swap(x1, x2);
                dy = y2 - y1;
                dx = x2 - x1;
            }
            m = dy / dx;
            float b = y1 - m*x1;
            for (float x = x1 + 1; x < x2; x++)
            {
                float y;
                y = round(x * m + b);
                pixel(x, y);
            }
        }
        else // Si Y > X, X por cada Y
        {
            if (y1 > y2)
            {
                swap(y1, y2);
                swap(x1, x2);
                dy = y2 - y1;
                dx = x2 - x1;
            }
            m = dy / dx;
            float b = y1 - m*x1;
            for (float y = y1 + 1; y < y2; y++)
            {
                float x;
                x = round((y-b) / m);
                pixel(x, y);
            }
        }
    }
    void linea_DDA(int x1, int y1, int x2, int y2) // Practica 2
    {
        pixel(x1, y1);
        if (x1 == x2 && y1 == x2) return;
        pixel(x2, y2);
        float dx, dy;
        dy = y2 - y1;
        dx = x2 - x1;
        if (abs((int)dx) >= abs((int)dy)) 
        {
            if (x1 > x2)
            {
                swap(y1, y2);
                swap(x1, x2);
                dy = y2 - y1;
                dx = x2 - x1;
            }
            float m = dy / dx;
            float y = y1;
            for (float x = x1 + 1; x < x2; x++)
            {
                y = y + m;
                pixel(x, round(y));
            }
        }
        else 
        {
            if (y1 > y2)
            {
                swap(y1, y2);
                swap(x1, x2);
                dy = y2 - y1;
                dx = x2 - x1;
            }
            float divm = dx / dy;
            float x = x1;
            for (float y = y1 + 1; y < y2; y++)
            {
                x = divm + x;
                pixel(round(x), y);
            }
        }
    }   
    // Esta es la linea de bresenham pero la modifique para poder interpolar dentro los x y originales
    vector<vector<double> > Linea_Bres(int x1, int y1, double ox1, double oy1, double z1, int x2, int y2, double ox2, double oy2, double z2)
    {
        vector<vector<double> > toReturn;
        vector<double> toInsert(5, 0);
        toInsert[0] = x1;
        toInsert[1] = y1;
        toInsert[2] = z1;
        toInsert[3] = ox1;
        toInsert[4] = oy1;
        toReturn.push_back(toInsert);
        if (x1 == x2 && y1 == x2)
            return toReturn;
        toInsert[0] = x2;
        toInsert[1] = y2;
        toInsert[2] = z2;
        toInsert[3] = ox2;
        toInsert[4] = oy2;
        toReturn.push_back(toInsert);

        int dx, dy;
        dy = y2 - y1;
        dx = x2 - x1;

        if (abs((int)dx) >= abs((int)dy))
        {
            if (x1 > x2) //Swap
            {
                swap(y1, y2);
                swap(x1, x2);
                swap(z1, z2);
                swap(ox1, ox2);
                swap(oy1, oy2);
                dy = y2 - y1;
                dx = x2 - x1;
            }
            int valor1 = 2 * abs((int)dy);
            int valor2 = valor1 - (2 * dx);

            int d = valor1 - dx;

            int y = y1;
            int add = (dy < 0) ? -1 : 1;

            double zstep = (z2 - z1) / abs(x2-x1);
            double zacum = z1 + zstep;

            double xstep = (ox2 - ox1) / abs(x2-x1);
            double xacum = ox1 + xstep;

            double ystep = (oy2 - oy1) / abs(x2-x1);
            double yacum = oy1 + ystep;

            for (int x = x1 + 1; x < x2; x++)
            {
                if (d < 0)
                    d += valor1;
                else
                {
                    d += valor2;
                    y += add;
                }
                toInsert[0] = x;
                toInsert[1] = y;
                toInsert[2] = zacum;
                toInsert[3] = xacum;
                toInsert[4] = yacum;
                toReturn.push_back(toInsert);
                zacum += zstep;
                xacum += xstep;
                yacum += ystep;
            }
        }
        else
        {
            if (y1 > y2) //Swap
            {
                swap(y1, y2);
                swap(x1, x2);
                swap(z1, z2);
                swap(ox1, ox2);
                swap(oy1, oy2);
                dy = y2 - y1;
                dx = x2 - x1;
            }

            int valor1 = 2 * abs((int)dx);
            int valor2 = valor1 - (2 * dy);

            int d = valor1 - dy;

            int x = x1;
            int add = (dx < 0) ? -1 : 1;

            double zstep = (z2 - z1) / abs(y2-y1);
            double zacum = z1 + zstep;

            double xstep = (ox2 - ox1) / abs(y2-y1);
            double xacum = ox1 + xstep;

            double ystep = (oy2 - oy1) / abs(y2-y1);
            double yacum = oy1 + ystep;

            for (int y = y1 + 1; y < y2; y++)
            {
                if (d < 0)
                    d += valor1;
                else
                {
                    d += valor2;
                    x += add;
                }
                toInsert[0] = x;
                toInsert[1] = y;
                toInsert[2] = zacum;
                toInsert[3] = xacum;
                toInsert[4] = yacum;
                toReturn.push_back(toInsert);
                zacum += zstep;
                xacum += xstep;
                yacum += ystep;
            }
        }
        return toReturn;
    }
    // Linea de bresenham pero ahora interpola los vertex normals para Phong
    // Phong calcula las normales en todos los puntos sobre las lineas
    vector<vector<double> > Linea_Bres_Phong(int x1, int y1, double ox1, double oy1, double z1, int x2, int y2, double ox2, double oy2, double z2, vector<double> n1, vector<double> n2)
    {
        vector<vector<double> > toReturn;
        vector<double> toInsert(8, 0);
        toInsert[0] = x1;
        toInsert[1] = y1;
        toInsert[2] = z1;
        toInsert[3] = ox1;
        toInsert[4] = oy1;
        toInsert[5] = n1[0];
        toInsert[6] = n1[1];
        toInsert[7] = n1[2];
        toReturn.push_back(toInsert);
        if (x1 == x2 && y1 == x2)
            return toReturn;
        toInsert[0] = x2;
        toInsert[1] = y2;
        toInsert[2] = z2;
        toInsert[3] = ox2;
        toInsert[4] = oy2;
        toInsert[5] = n2[0];
        toInsert[6] = n2[1];
        toInsert[7] = n2[2];
        toReturn.push_back(toInsert);

        int dx, dy;
        dy = y2 - y1;
        dx = x2 - x1;

        if (abs((int)dx) >= abs((int)dy))
        {
            if (x1 > x2) //Swap
            {
                swap(y1, y2);
                swap(x1, x2);
                swap(z1, z2);
                swap(ox1, ox2);
                swap(oy1, oy2);
                swap(n1, n2);
                dy = y2 - y1;
                dx = x2 - x1;
            }
            int valor1 = 2 * abs((int)dy);
            int valor2 = valor1 - (2 * dx);

            int d = valor1 - dx;

            int y = y1;
            int add = (dy < 0) ? -1 : 1;

            double zstep = (z2 - z1) / abs(x2-x1);
            double zacum = z1 + zstep;

            double xstep = (ox2 - ox1) / abs(x2-x1);
            double xacum = ox1 + xstep;

            double ystep = (oy2 - oy1) / abs(x2-x1);
            double yacum = oy1 + ystep;

            double xnstep = (n2[0] - n1[0]) / abs(x2 - x1);
            double xnacum = n1[0] + xnstep;
            double ynstep = (n2[1] - n1[1]) / abs(x2 - x1);
            double ynacum = n1[1] + ynstep;
            double znstep = (n2[2] - n1[2]) / abs(x2 - x1);
            double znacum = n1[2] + znstep;

            
            for (int x = x1 + 1; x < x2; x++)
            {
                if (d < 0)
                    d += valor1;
                else
                {
                    d += valor2;
                    y += add;
                }
                toInsert[0] = x;
                toInsert[1] = y;
                toInsert[2] = zacum;
                toInsert[3] = xacum;
                toInsert[4] = yacum;
                toInsert[5] = xnacum;
                toInsert[6] = ynacum;
                toInsert[7] = znacum;
                toReturn.push_back(toInsert);
                zacum += zstep;
                xacum += xstep;
                yacum += ystep;
                xnacum += xnstep;
                ynacum += ynstep;
                znacum += znstep;
            }
        }
        else
        {
            if (y1 > y2) //Swap
            {
                swap(y1, y2);
                swap(x1, x2);
                swap(z1, z2);
                swap(ox1, ox2);
                swap(oy1, oy2);
                swap(n1, n2);
                dy = y2 - y1;
                dx = x2 - x1;
            }

            int valor1 = 2 * abs((int)dx);
            int valor2 = valor1 - (2 * dy);

            int d = valor1 - dy;

            int x = x1;
            int add = (dx < 0) ? -1 : 1;

            double zstep = (z2 - z1) / abs(y2-y1);
            double zacum = z1 + zstep;

            double xstep = (ox2 - ox1) / abs(y2-y1);
            double xacum = ox1 + xstep;

            double ystep = (oy2 - oy1) / abs(y2-y1);
            double yacum = oy1 + ystep;

            double xnstep = (n2[0] - n1[0]) / abs(y2 - y1);
            double xnacum = n1[0] + xnstep;
            double ynstep = (n2[1] - n1[1]) / abs(y2 - y1);
            double ynacum = n1[1] + ynstep;
            double znstep = (n2[2] - n1[2]) / abs(y2 - y1);
            double znacum = n1[2] + znstep;

            for (int y = y1 + 1; y < y2; y++)
            {
                if (d < 0)
                    d += valor1;
                else
                {
                    d += valor2;
                    x += add;
                }
                toInsert[0] = x;
                toInsert[1] = y;
                toInsert[2] = zacum;
                toInsert[3] = xacum;
                toInsert[4] = yacum;
                toInsert[5] = xnacum;
                toInsert[6] = ynacum;
                toInsert[7] = znacum;
                toReturn.push_back(toInsert);
                zacum += zstep;
                xacum += xstep;
                yacum += ystep;
                xnacum += xnstep;
                ynacum += ynstep;
                znacum += znstep;
            }
        }
        return toReturn;
    }
    // Linea de bresenham pero ahora interpola colores para Gouraud
    // Gouraud calcula los colores finales de los vertices y los interpola sobre los lados
    vector<vector<double> > Linea_Bres_Gouraud(int x1, int y1, double z1, int x2, int y2, double z2, vector<double> n1, vector<double> n2)
    {
        vector<vector<double> > toReturn;
        vector<double> toInsert(6, 0);
        toInsert[0] = x1;
        toInsert[1] = y1;
        toInsert[2] = z1;
        toInsert[3] = n1[0];
        toInsert[4] = n1[1];
        toInsert[5] = n1[2];
        toReturn.push_back(toInsert);
        if (x1 == x2 && y1 == x2)
            return toReturn;
        toInsert[0] = x2;
        toInsert[1] = y2;
        toInsert[2] = z2;
        toInsert[3] = n2[0];
        toInsert[4] = n2[1];
        toInsert[5] = n2[2];
        toReturn.push_back(toInsert);

        int dx, dy;
        dy = y2 - y1;
        dx = x2 - x1;

        if (abs((int)dx) >= abs((int)dy))
        {
            if (x1 > x2) //Swap
            {
                swap(y1, y2);
                swap(x1, x2);
                swap(z1, z2);
                swap(n1, n2);
                dy = y2 - y1;
                dx = x2 - x1;
            }
            int valor1 = 2 * abs((int)dy);
            int valor2 = valor1 - (2 * dx);

            int d = valor1 - dx;

            int y = y1;
            int add = (dy < 0) ? -1 : 1;

            double zstep = (z2 - z1) / abs(x2 - x1);
            double zacum = z1 + zstep;

            double xnstep = (n2[0] - n1[0]) / abs(x2 - x1);
            double xnacum = n1[0] + xnstep; //R
            double ynstep = (n2[1] - n1[1]) / abs(x2 - x1);
            double ynacum = n1[1] + ynstep; //G
            double znstep = (n2[2] - n1[2]) / abs(x2 - x1);
            double znacum = n1[2] + znstep; //B

            for (int x = x1 + 1; x < x2; x++)
            {
                if (d < 0)
                    d += valor1;
                else
                {
                    d += valor2;
                    y += add;
                }
                toInsert[0] = x;
                toInsert[1] = y;
                toInsert[2] = zacum;
                toInsert[3] = xnacum;
                toInsert[4] = ynacum;
                toInsert[5] = znacum;
                toReturn.push_back(toInsert);
                zacum += zstep;
                xnacum += xnstep;
                ynacum += ynstep;
                znacum += znstep;
            }
        }
        else
        {
            if (y1 > y2) //Swap
            {
                swap(y1, y2);
                swap(x1, x2);
                swap(z1, z2);
                swap(n1, n2);
                dy = y2 - y1;
                dx = x2 - x1;
            }

            int valor1 = 2 * abs((int)dx);
            int valor2 = valor1 - (2 * dy);

            int d = valor1 - dy;

            int x = x1;
            int add = (dx < 0) ? -1 : 1;

            double zstep = (z2 - z1) / abs(y2 - y1);
            double zacum = z1 + zstep;

            double xnstep = (n2[0] - n1[0]) / abs(y2 - y1);
            double xnacum = n1[0] + xnstep; //R
            double ynstep = (n2[1] - n1[1]) / abs(y2 - y1);
            double ynacum = n1[1] + ynstep; //G
            double znstep = (n2[2] - n1[2]) / abs(y2 - y1);
            double znacum = n1[2] + znstep; //B

            for (int y = y1 + 1; y < y2; y++)
            {
                if (d < 0)
                    d += valor1;
                else
                {
                    d += valor2;
                    x += add;
                }
                toInsert[0] = x;
                toInsert[1] = y;
                toInsert[2] = zacum;
                toInsert[3] = xnacum;
                toInsert[4] = ynacum;
                toInsert[5] = znacum;
                toReturn.push_back(toInsert);
                zacum += zstep;
                xnacum += xnstep;
                ynacum += ynstep;
                znacum += znstep;
            }
        }
        return toReturn;
    }
    
    void iluminar_Luz(Luz luz)
    {
        vector<double> vectorZ1(3, 0);
        vectorZ1[2] = 1; // Vector z para la camara
        for (int iy = 0; iy < ry; iy++)
        {
            for (int ix = 0; ix < rx; ix++)
            {
                if (xyPixel[ix][iy][0] != (double)INFINITY && xyPixel[ix][iy][1] != (double)INFINITY)
                {
                    Material mt = materialPixel[ix][iy];
                    vector<double> pnormal = normalPixel[ix][iy];
                    vector<double> val = mt.obtenerValores();
                    double ak, dk, sk, ns; // atenuación, difuso y especular
                    ak = val[0];
                    dk = val[1];
                    sk = val[2];
                    ns = val[3];

                    vector<double> lightPosition, attenVar, pixpos; //Posicion de luz, variable de atenuación, posicion original del pixel
                    vector<int> lightColor;
                    pixpos = vector<double>(3, 0);
                    pixpos[0] = xyPixel[ix][iy][0];
                    pixpos[1] = xyPixel[ix][iy][1];
                    pixpos[2] = zBuffer[ix][iy];
                    vector<int> newColor = vector<int>(3, 0);
                    newColor[0] = image[ix][iy][0];
                    newColor[1] = image[ix][iy][1];
                    newColor[2] = image[ix][iy][2];

                    lightPosition = luz.regresarPosicion();
                    lightColor = luz.regresarColor();
                    attenVar = luz.regresarAtenuacion();
                    double angle = 0; 
                    double coneAngle = 10;
                    vector<double> dir;
                    if (luz.is_Cone()) // Si la luz es un cono
                    {
                        vector<double> paraAngulo = op.restarVector(pixpos, lightPosition); //  Se compara el ángulo entre el origen y el punto del pixel con él ángulo de la luz
                        dir = luz.regresarDireccion();
                        angle = acos((op.productoPunto(paraAngulo, dir)) / (op.magnitud_Vector(paraAngulo) * op.magnitud_Vector(dir))) * (180 / M_PI);
                        coneAngle = luz.regresarAngulo();
                    }
                    else
                    {
                        dir = op.restarVector(pixpos, lightPosition);
                    }
                    if (angle < coneAngle) // Si el ángulo es menor al ángulo de apertura, se ilumina
                    {
                        vector<double> reflejo = op.reflejoVector(dir, pnormal);
                        double diffAngle = acos((op.productoPunto(pnormal, reflejo)) / (op.magnitud_Vector(pnormal) * op.magnitud_Vector(reflejo))); // Luz Difusa
                        double diffValue = dk * cos(diffAngle);
                        if(diffValue < 0 ) 
                            diffValue = 0;
                        double specAngle = acos((op.productoPunto(vectorZ1, reflejo)) / (op.magnitud_Vector(vectorZ1) * op.magnitud_Vector(reflejo))); // Luz especular
                        double specValue = sk * pow(cos(specAngle), ns);

                        double d = op.magnitud_Vector(op.restarVector(pixpos, lightPosition));
                        double attenuation = 1 / ((attenVar[0] * d * d) + (attenVar[1] * d) + attenVar[2]);
                        for(int i=0 ; i < 3; i++)
                        {
                            newColor[i] += lightColor[i] * attenuation * diffValue;
                            newColor[i] += lightColor[i] * attenuation * specValue;
                        }
                    }
                    newColor[0] += ak * luzAmbiental[0]; // Luz ambiental
                    newColor[1] += ak * luzAmbiental[1]; 
                    newColor[2] += ak * luzAmbiental[2]; 
                    for (int cnt = 0; cnt < 3; cnt++)
                    {
                        if (newColor[cnt] > 255)
                            newColor[cnt] = 255;
                        if (newColor[cnt] < 0)
                            newColor[cnt] = 0;
                    }
                    image[ix][iy][0] = newColor[0];
                    image[ix][iy][1] = newColor[1];
                    image[ix][iy][2] = newColor[2];
                }
            }
        }
    }
    
    void pintarOBJPlano(OBJ objeto, Luz luz)
    {
        OBJ obj = objeto;
        vector<double> vectorZ1(3,0); //Vector de la camara
        vector<vector<double> > normalesObj = obj.regresarNormales();
        // Proyectar
        Transformacion proyeccion(0);
        vector<vector<double> > tabla(4, vector<double>(4, 0));
        for (int i = 0; i < 3; i++)
           tabla[i][i] = 1;
        tabla[3][2] = 1 / f;
        proyeccion.definir_Tabla(tabla);
        obj.aplicar_Transformacion(proyeccion);
        obj.dividir_Homogeneo();
        vectorZ1[2] = 1; // Vector z para la camara 
        
        vector<vector<double> > puntosObj = obj.regresarPuntos();
        vector<vector<double> > ladosObj = obj.regresarLados();
        vector<vector<double> > carasObj = obj.regresarCaras();
        vector<vector<double> > puntosOriginales = objeto.regresarPuntos();
        Material materialObj = obj.regresarMaterial();

        for (int c = 0; c < carasObj.size(); c++) // Por cada cara
        {
            // Se calcula el ángulo entre la cámara ya la cara
            double angle = acos((op.productoPunto(normalesObj[c], vectorZ1))/(sqrt(pow(normalesObj[c][0], 2) + pow(normalesObj[c][1], 2) + pow(normalesObj[c][2], 2)))) * (180 / M_PI); //La camara estan normalizados.
            if(180 >= angle && angle >= 90) 
            {
                // Pintamos la cara
                vector<vector<double> > toPaint;
                vector<vector<double> > original;
                for (int l = 0; l < 3; l++) // Se obtienen los vertices
                {
                    int x1, y1, x2, y2;
                    double z1, z2;
                    double ox1, ox2, oy1, oy2; // Se sacan los puntos del lado
                    x1 = round(puntosObj[ladosObj[carasObj[c][l] - 1][0] - 1][0]);
                    y1 = round(puntosObj[ladosObj[carasObj[c][l] - 1][0] - 1][1]);
                    z1 = puntosOriginales[ladosObj[carasObj[c][l] - 1][0] - 1][2];
                    x2 = round(puntosObj[ladosObj[carasObj[c][l] - 1][1] - 1][0]);
                    y2 = round(puntosObj[ladosObj[carasObj[c][l] - 1][1] - 1][1]);
                    z2 = puntosOriginales[ladosObj[carasObj[c][l] - 1][1] - 1][2];

                    ox1 = puntosOriginales[ladosObj[carasObj[c][l] - 1][0] - 1][0];
                    oy1 = puntosOriginales[ladosObj[carasObj[c][l] - 1][0] - 1][1];
                    ox2 = puntosOriginales[ladosObj[carasObj[c][l] - 1][1] - 1][0];
                    oy2 = puntosOriginales[ladosObj[carasObj[c][l] - 1][1] - 1][1];
                    // Se obtienen los puntos de la línea
                    vector<vector<double> > puntos = Linea_Bres(x1, y1, ox1, oy1, z1, x2, y2, ox2, oy2, z2);
                    for(vector<double> p : puntos)
                        toPaint.push_back(p);
                }
                // Ordenamiento de burbuja para ordenar las Y
                int i, j;
                for (i = 0; i < toPaint.size() - 1; i++)
                {
                    for (j = 0; j < toPaint.size() - i - 1; j++)
                    {
                        if (toPaint[j][1] > toPaint[j + 1][1])
                        {
                            vector<double> aux = toPaint[j];
                            toPaint[j] = toPaint[j + 1];
                            toPaint[j + 1] = aux;
                        }
                    }
                }
                int yRev = toPaint[0][1];
                vector<vector<double> > bulk; //Guarda los puntos que comparten la misma y
                for (int i = 0; i < toPaint.size()+1; i++)
                {
                    if( i < toPaint.size() && yRev == toPaint[i][1]) { bulk.push_back(toPaint[i]); }
                    else
                    { // Una vez encuentre un punto con una y diferente, comienza a rellenar sobre esa y
                        int xmin = rx; // saca la x máxima y minima de la línea que se pintara
                        int xmax = -1;
                        double zmin, zmax, oxmin, oxmax, oymin, oymax;
                        for (int i2 = 0; i2 < bulk.size(); i2++)
                        {
                            if(xmin > bulk[i2][0]) 
                            { 
                                xmin = bulk[i2][0]; 
                                zmin = bulk[i2][2];
                                oxmin = bulk[i2][3];
                                oymin = bulk[i2][4];
                            }
                            if(xmax < bulk[i2][0]) 
                            { 
                                xmax = bulk[i2][0]; 
                                zmax = bulk[i2][2];
                                oxmax = bulk[i2][3];
                                oymax = bulk[i2][4];
                            }
                        }
                        pixel(xmin, yRev, oxmin, oymin, zmin, obj.regresarRGB(), materialObj, normalesObj[c]); // Esquina izquierda
                        if(xmin != xmax)
                        {
                            pixel(xmax, yRev, oxmax, oymax, zmax, obj.regresarRGB(), materialObj, normalesObj[c]);  // Esquina derecha

                            double stepz = (zmax - zmin) / abs(xmax - xmin);
                            double acumuladoz = zmin + stepz;

                            double xstep = (oxmax - oxmin) / abs(xmax - xmin);
                            double xacum = oxmin + xstep;

                            double ystep = (oymax - oymin) / abs(xmax - xmin);
                            double yacum = oymin + ystep;

                            for (int cnt = xmin+1; cnt < xmax; cnt++)
                            {
                                pixel(cnt, yRev, xacum, yacum, acumuladoz, obj.regresarRGB(), materialObj, normalesObj[c]);
                                acumuladoz += stepz;
                                xacum += xstep;
                                yacum += ystep;
                            }
                        }
                        if( i < toPaint.size()){ // Cambiar la Y a revisar
                            yRev = toPaint[i][1];
                            bulk.clear();
                            bulk.push_back(toPaint[i]);
                        }
                    }
                }
            }
        }
        iluminar_Luz(luz); // Vamos a iluminar sobre de los pixeles
        return;
    }  
    void pintarOBJPhong(OBJ objeto, Luz luz)
    {
        OBJ obj = objeto;
        vector<double> vectorZ1(3,0); //Vector de la camara
        vector<vector<double> > normalesObj = obj.regresarNormales();
        
        // Proyectar
        Transformacion proyeccion(0);
        vector<vector<double> > tabla(4, vector<double>(4, 0));
        for (int i = 0; i < 3; i++)
            tabla[i][i] = 1;
        tabla[3][2] = 1 / f;
        proyeccion.definir_Tabla(tabla);
        obj.aplicar_Transformacion(proyeccion);
        obj.dividir_Homogeneo();

        vectorZ1[2] = 1; // Vector z para la camara
        
        vector<vector<double> > puntosObj = obj.regresarPuntos();
        vector<vector<double> > ladosObj = obj.regresarLados();
        vector<vector<double> > carasObj = obj.regresarCaras();
        vector<vector<double> > puntosOriginales = objeto.regresarPuntos();
        Material materialObj = obj.regresarMaterial();
        //Calculate vertex Normals
        vector<vector<double> > vertexNormal = objeto.regresarNormalesVertices();
        for (int c = 0; c < carasObj.size(); c++)
        {
            double angle = acos((op.productoPunto(normalesObj[c], vectorZ1))/(sqrt(pow(normalesObj[c][0], 2) + pow(normalesObj[c][1], 2) + pow(normalesObj[c][2], 2)))) * (180 / M_PI); //La camara estan normalizados.
            //checar angulo
            if(180 >= angle && angle >= 90) //Show it! 
            {
                vector<vector<double> > toPaint;
                for (int l = 0; l < 3; l++)
                {
                    int x1, y1, x2, y2;
                    double z1, z2;
                    double ox1, ox2, oy1, oy2; 
                    x1 = round(puntosObj[ladosObj[carasObj[c][l] - 1][0] - 1][0]);
                    y1 = round(puntosObj[ladosObj[carasObj[c][l] - 1][0] - 1][1]);
                    z1 = puntosOriginales[ladosObj[carasObj[c][l] - 1][0] - 1][2];
                    x2 = round(puntosObj[ladosObj[carasObj[c][l] - 1][1] - 1][0]);
                    y2 = round(puntosObj[ladosObj[carasObj[c][l] - 1][1] - 1][1]);
                    z2 = puntosOriginales[ladosObj[carasObj[c][l] - 1][1] - 1][2];

                    ox1 = puntosOriginales[ladosObj[carasObj[c][l] - 1][0] - 1][0];
                    oy1 = puntosOriginales[ladosObj[carasObj[c][l] - 1][0] - 1][1];
                    ox2 = puntosOriginales[ladosObj[carasObj[c][l] - 1][1] - 1][0];
                    oy2 = puntosOriginales[ladosObj[carasObj[c][l] - 1][1] - 1][1];

                    vector<double> vn1 = vertexNormal[ladosObj[carasObj[c][l] - 1][0] - 1];
                    vector<double> vn2 = vertexNormal[ladosObj[carasObj[c][l] - 1][1] - 1];

                    vector<vector<double> > puntos = Linea_Bres_Phong(x1, y1, ox1, oy1, z1, x2, y2, ox2, oy2, z2, vn1, vn2);
                    // [3][4] = xy originales  ||   [5][6][7] = normal del pixel
                    for(vector<double> p : puntos)
                        toPaint.push_back(p);
                }
                // Ordenamiento de burbuja
                int i, j;
                for (i = 0; i < toPaint.size() - 1; i++)
                {
                    for (j = 0; j < toPaint.size() - i - 1; j++)
                    {
                        if (toPaint[j][1] > toPaint[j + 1][1])
                        {
                            vector<double> aux = toPaint[j];
                            toPaint[j] = toPaint[j + 1];
                            toPaint[j + 1] = aux;
                        }
                    }
                }
                int yRev = toPaint[0][1];
                vector<vector<double> > bulk;
                for (int i = 0; i < toPaint.size()+1; i++)
                {
                    if( i < toPaint.size() && yRev == toPaint[i][1]) { bulk.push_back(toPaint[i]); }
                    else
                    {
                         
                        int xmin = rx;
                        int xmax = -1;
                        double zmin, zmax, oxmin, oxmax, oymin, oymax;
                        double xnmin, xnmax, ynmin, ynmax, znmin, znmax; // Interpolando las normales
                        for (int i2 = 0; i2 < bulk.size(); i2++)
                        {
                            if(xmin > bulk[i2][0]) // A diferencia del plano, aquí se interpolan las nomrales
                            { 
                                xmin = bulk[i2][0]; 
                                zmin = bulk[i2][2];
                                oxmin = bulk[i2][3];
                                oymin = bulk[i2][4];
                                xnmin = bulk[i2][5];
                                ynmin = bulk[i2][6];
                                znmin = bulk[i2][7];
                            }
                            if(xmax < bulk[i2][0]) 
                            { 
                                xmax = bulk[i2][0]; 
                                zmax = bulk[i2][2];
                                oxmax = bulk[i2][3];
                                oymax = bulk[i2][4];
                                xnmax = bulk[i2][5];
                                ynmax = bulk[i2][6];
                                znmax = bulk[i2][7];
                            }
                        }
                        pixel(xmin, yRev, oxmin, oymin, zmin, obj.regresarRGB(), materialObj, {xnmin, ynmin, znmin});
                        if(xmin != xmax)
                        {
                            pixel(xmax, yRev, oxmax, oymax, zmax, obj.regresarRGB(), materialObj, {xnmax, ynmax, znmax});

                            double stepz = (zmax - zmin) / abs(xmax - xmin);
                            double acumuladoz = zmin + stepz;
                            double xstep = (oxmax - oxmin) / abs(xmax - xmin);
                            double xacum = oxmin + xstep;
                            double ystep = (oymax - oymin) / abs(xmax - xmin);
                            double yacum = oymin + ystep;

                            double xnstep = (xnmax - xnmin) / abs(xmax - xmin);
                            double xnacum = xnmin + xnstep;
                            double ynstep = (ynmax - ynmin) / abs(xmax - xmin);
                            double ynacum = ynmin + ynstep;
                            double znstep = (znmax - znmin) / abs(xmax - xmin);
                            double znacum = znmin + znstep;

                            for (int cnt = xmin+1; cnt < xmax; cnt++)
                            {
                                pixel(cnt, yRev, xacum, yacum, acumuladoz, obj.regresarRGB(), materialObj, {xnacum, ynacum, znacum});
                                acumuladoz += stepz;
                                xacum += xstep;
                                yacum += ystep;
                                xnacum += xnstep;
                                ynacum += ynstep;
                                znacum += znstep;
                            }
                        }
                        if(i < toPaint.size()){
                            //Change yrev 
                            yRev = toPaint[i][1];
                            bulk.clear();
                            bulk.push_back(toPaint[i]);
                        }
                    }
                }
            }
        }
        iluminar_Luz(luz); // Vamos a iluminar encima de los pixeles
        return;
    }
    void pintarOBJGouraud(OBJ objeto, Luz luz)
    {
        vector<double> vectorZ1(3, 0);        //Vector de la camara
         //Suponiendo 'z'
        vectorZ1[2] = 1; // Vector z para la camara

        OBJ obj = objeto;
        vector<vector<double> > normalesObj = obj.regresarNormales();
        Transformacion proyeccion(0);
        vector<vector<double> > tabla(4, vector<double>(4, 0));
        for (int i = 0; i < 3; i++)
            tabla[i][i] = 1;
        tabla[3][2] = 1 / f;
        proyeccion.definir_Tabla(tabla);
        obj.aplicar_Transformacion(proyeccion);
        obj.dividir_Homogeneo();

        vector<vector<double> > puntosObj = obj.regresarPuntos();
        vector<vector<double> > ladosObj = obj.regresarLados();
        vector<vector<double> > carasObj = obj.regresarCaras();
        vector<vector<double> > puntosOriginales = objeto.regresarPuntos();
        Material materialObj = obj.regresarMaterial();
        vector<vector<double> > vertexNormal = objeto.regresarNormalesVertices(); //Calcula las normales de los vertices
        for (int c = 0; c < carasObj.size(); c++)
        {
            double angle = acos((op.productoPunto(normalesObj[c], vectorZ1)) / (sqrt(pow(normalesObj[c][0], 2) + pow(normalesObj[c][1], 2) + pow(normalesObj[c][2], 2)))) * (180 / M_PI); //La camara estan normalizados.
            //Si se encuentra dentro del ángulo de apertura del cono, lo muestra
            if (180 >= angle && angle >= 90) 
            {
                vector<vector<double> > toPaint;
                for (int l = 0; l < 3; l++)
                {
                    int x1, y1, x2, y2;
                    double z1, z2;
                    double ox1, ox2, oy1, oy2; //o stands for original (before projecting)
                    x1 = round(puntosObj[ladosObj[carasObj[c][l] - 1][0] - 1][0]);
                    y1 = round(puntosObj[ladosObj[carasObj[c][l] - 1][0] - 1][1]);
                    z1 = puntosOriginales[ladosObj[carasObj[c][l] - 1][0] - 1][2];
                    x2 = round(puntosObj[ladosObj[carasObj[c][l] - 1][1] - 1][0]);
                    y2 = round(puntosObj[ladosObj[carasObj[c][l] - 1][1] - 1][1]);
                    z2 = puntosOriginales[ladosObj[carasObj[c][l] - 1][1] - 1][2];

                    ox1 = puntosOriginales[ladosObj[carasObj[c][l] - 1][0] - 1][0];
                    oy1 = puntosOriginales[ladosObj[carasObj[c][l] - 1][0] - 1][1];
                    ox2 = puntosOriginales[ladosObj[carasObj[c][l] - 1][1] - 1][0];
                    oy2 = puntosOriginales[ladosObj[carasObj[c][l] - 1][1] - 1][1];

                    vector<double> vn1 = vertexNormal[ladosObj[carasObj[c][l] - 1][0] - 1];
                    vector<double> vn2 = vertexNormal[ladosObj[carasObj[c][l] - 1][1] - 1];
                    Material mt = obj.regresarMaterial();

                    vector<double> val = mt.obtenerValores();
                    double ak, dk, sk, ns;
                    ak = val[0];
                    dk = val[1];
                    sk = val[2];
                    ns = val[3];

                    vector<double> lightPosition, attenVar, pixpos, pixpos2;
                    vector<int> lightColor;
                    pixpos = vector<double>(3, 0);
                    pixpos2 = vector<double>(3, 0);
                    pixpos[0] = ox1;
                    pixpos[1] = oy1;
                    pixpos[2] = z1;
                    pixpos2[0] = ox2;
                    pixpos2[1] = oy2;
                    pixpos2[2] = z2;
                    vector<double> newColorV1 = vector<double>(3, 0);
                    vector<double> newColorV2 = vector<double>(3, 0);
                    vector<int> colorOriginal = obj.regresarRGB();
                    for (int d = 0; d < 3; d++)
                    {
                        newColorV1[d] = (colorOriginal[d]);
                        newColorV2[d] = (colorOriginal[d]);
                    }

                    lightPosition = luz.regresarPosicion();
                    lightColor = luz.regresarColor();
                    attenVar = luz.regresarAtenuacion();
                    double angle = 0;
                    double coneAngle = 10;
                    vector<double> dir;
                    bool v1_in = false, v2_in = false;
                    if (luz.is_Cone())
                    {
                        vector<double> paraAngulo = op.restarVector(pixpos, lightPosition);
                        dir = luz.regresarDireccion();
                        angle = acos((op.productoPunto(paraAngulo, dir)) / (op.magnitud_Vector(paraAngulo) * op.magnitud_Vector(dir))) * (180 / M_PI);
                        coneAngle = luz.regresarAngulo();
                        if (angle < coneAngle)
                            v1_in = true;
                        paraAngulo = op.restarVector(pixpos2, lightPosition);
                        angle = acos((op.productoPunto(paraAngulo, dir)) / (op.magnitud_Vector(paraAngulo) * op.magnitud_Vector(dir))) * (180 / M_PI);
                        if (angle < coneAngle)
                            v2_in = true;
                    }
                    else
                    {
                        dir = op.restarVector(pixpos, lightPosition);
                        v1_in = true;
                        v2_in = true;
                    }
                    // Pain vertex 1
                    if (v1_in)
                    {
                        vector<double> reflejo = op.reflejoVector(dir, vn1);
                        double diffAngle = acos((op.productoPunto(vn1, reflejo)) / (op.magnitud_Vector(vn1) * op.magnitud_Vector(reflejo)));
                        double diffValue = dk * cos(diffAngle);
                        if (diffValue < 0)
                            diffValue = 0;
                        double specAngle = acos((op.productoPunto(vectorZ1, reflejo)) / (op.magnitud_Vector(vectorZ1) * op.magnitud_Vector(reflejo)));
                        double specValue = sk * pow(cos(specAngle), ns);
                        double d = op.magnitud_Vector(op.restarVector(pixpos, lightPosition));
                        double attenuation = 1 / ((attenVar[0] * d * d) + (attenVar[1] * d) + attenVar[2]);
                        for (int i = 0; i < 3; i++)
                        {
                            newColorV1[i] += lightColor[i] * attenuation * diffValue;
                            newColorV1[i] += lightColor[i] * attenuation * specValue;
                        }
                    }
                    // Pintar el segundo vertice
                    if (v2_in)
                    {
                        vector<double> reflejo = op.reflejoVector(dir, vn2);
                        double diffAngle = acos((op.productoPunto(vn2, reflejo)) / (op.magnitud_Vector(vn2) * op.magnitud_Vector(reflejo)));
                        double diffValue = dk * cos(diffAngle);
                        if (diffValue < 0)
                            diffValue = 0;
                        double specAngle = acos((op.productoPunto(vectorZ1, reflejo)) / (op.magnitud_Vector(vectorZ1) * op.magnitud_Vector(reflejo)));
                        double specValue = sk * pow(cos(specAngle), ns);
                        double d = op.magnitud_Vector(op.restarVector(pixpos2, lightPosition));
                        double attenuation = 1 / ((attenVar[0] * d * d) + (attenVar[1] * d) + attenVar[2]);
                        for (int i = 0; i < 3; i++)
                        {
                            newColorV2[i] += lightColor[i] * attenuation * diffValue;
                            newColorV2[i] += lightColor[i] * attenuation * specValue;
                        }
                    }
                    newColorV1[0] += ak * luzAmbiental[0]; // Luz ambiental en el vertice 1
                    newColorV1[1] += ak * luzAmbiental[1]; 
                    newColorV1[2] += ak * luzAmbiental[2]; 
                    newColorV2[0] += ak * luzAmbiental[0]; // Luz ambiental en el vertice 2
                    newColorV2[1] += ak * luzAmbiental[1]; 
                    newColorV2[2] += ak * luzAmbiental[2]; 
                    for (int cnt = 0; cnt < 3; cnt++)
                    {
                        if (newColorV1[cnt] > 255)
                            newColorV1[cnt] = 255;
                        if (newColorV2[cnt] > 255)
                            newColorV2[cnt] = 255;
                    }
                    vector<vector<double> > puntos = Linea_Bres_Gouraud(x1, y1, z1, x2, y2, z2, newColorV1, newColorV2);
                    // [4][5][6]  = color del pixel interpolado
                    for (vector<double> p : puntos)
                        toPaint.push_back(p);
                }
                // Ordenamiento de burbuja
                int i, j;
                for (i = 0; i < toPaint.size() - 1; i++)
                {
                    for (j = 0; j < toPaint.size() - i - 1; j++)
                    {
                        if (toPaint[j][1] > toPaint[j + 1][1])
                        {
                            vector<double> aux = toPaint[j];
                            toPaint[j] = toPaint[j + 1];
                            toPaint[j + 1] = aux;
                        }
                    }
                }
                int yRev = toPaint[0][1];
                vector<vector<double> > bulk;
                for (int i = 0; i < toPaint.size() + 1; i++)
                {
                    if ( i < toPaint.size() && yRev == toPaint[i][1])
                    {
                        bulk.push_back(toPaint[i]);
                    }
                    else
                    {
                         
                        int xmin = rx;
                        int xmax = -1;
                        double zmin, zmax;
                        double xnmin, xnmax, ynmin, ynmax, znmin, znmax; // Interpolando los colores
                        for (int i2 = 0; i2 < bulk.size(); i2++)
                        {
                            if (xmin > bulk[i2][0])
                            {
                                xmin = bulk[i2][0];
                                zmin = bulk[i2][2];
                                xnmin = bulk[i2][3];
                                ynmin = bulk[i2][4];
                                znmin = bulk[i2][5];
                            }
                            if (xmax < bulk[i2][0])
                            {
                                xmax = bulk[i2][0];
                                zmax = bulk[i2][2];
                                xnmax = bulk[i2][3];
                                ynmax = bulk[i2][4];
                                znmax = bulk[i2][5];
                            }
                        }
                        pixel(xmin, yRev, zmin, {(int)xnmin, (int)ynmin, (int)znmin});
                        if (xmin != xmax)
                        {
                            pixel(xmax, yRev, zmax, {(int)xnmax, (int)ynmax, (int)znmax});

                            double stepz = (zmax - zmin) / abs(xmax - xmin);
                            double acumuladoz = zmin + stepz;

                            double xnstep = (xnmax - xnmin) / abs(xmax - xmin);
                            double xnacum = xnmin + xnstep; //R
                            double ynstep = (ynmax - ynmin) / abs(xmax - xmin);
                            double ynacum = ynmin + ynstep; //G
                            double znstep = (znmax - znmin) / abs(xmax - xmin);
                            double znacum = znmin + znstep; //B

                            for (int cnt = xmin + 1; cnt < xmax; cnt++)
                            {
                                pixel(cnt, yRev, acumuladoz, {(int)xnacum, (int)ynacum, (int)znacum});
                                acumuladoz += stepz;
                                xnacum += xnstep;
                                ynacum += ynstep;
                                znacum += znstep;
                            }
                        }
                        if (i < toPaint.size())
                        {
                            //Change yrev
                            yRev = toPaint[i][1];
                            bulk.clear();
                            bulk.push_back(toPaint[i]);
                        }
                    }
                }
            }
        }
        return;
    }

    void setLuzAmbiental(int _r, int _g, int _b) { 
        luzAmbiental[0] = _r;
        luzAmbiental[1] = _g;
        luzAmbiental[2] = _b;
    }
};

int main(int argc, const char** argv) {

    // Vamos a probar el código 

    // Creamos el raster
    Raster imagen(1920,1080); // Tamano Full HD
    imagen.setTamano_Focal(8000); // Le Asignamos un nuevo tamano focal a la camara
    imagen.setLuzAmbiental(5, 5, 5); // Le Asignamos un valor a la luz ambiental

    OBJ monkey("Monkey", false); // Creamos al objeto monkey, lo acomodeamos y lo convertimos en un vlf
    Transformacion rotar('y', M_PI);
    Transformacion mover(2300,1350,20000);
    Transformacion escalar(200);
    monkey.aplicar_Transformacion(rotar);
    monkey.aplicar_Transformacion(escalar);
    monkey.aplicar_Transformacion(mover); 
    monkey.transformar_VLF();
    monkey.cambiar_color({255,0,0}); // Le ponemos un color rojizo
    //monkey.cambiar_Material(Material(1,0.7,0.01,100));

    Luz luz(                    // Creamos la luz que iluminara la escena
        {-500, -100, 13000}, 
        {255, 255, 255},
        {0.00000001,0.00001,0}
    );

    // Vamos a copiar el raster para crear una imagen plana , una imagen de Phong y una de Gouraud
    Raster imagen_plana, imagen_phong, imagen_gouraud;
    
    // Para la imagen plana
    imagen_plana = imagen;
    imagen_plana.pintarOBJPlano(monkey, luz);
    imagen_plana.crearPPM("Monkey_Plano");
    cout << "Completada imagen plana :)";
    // Para la imagen de phong
    imagen_phong = imagen;
    imagen_phong.pintarOBJPhong(monkey, luz);
    imagen_phong.crearPPM("Monkey_Phong");
    cout << "\nCompletada imagen phong :)";
    // Para la imagen de gouraud
    imagen_gouraud = imagen;
    imagen_gouraud.pintarOBJGouraud(monkey, luz);
    imagen_gouraud.crearPPM("Monkey_Gouraud");
    cout << "\nCompletada imagen gouraud :)\n";

    return 0;
}
