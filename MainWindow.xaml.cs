using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Xml.Linq;

namespace Electric
{
    /// <summary>
    /// Логика взаимодействия для MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        double a = 0; // Only for drawning
        double kr; // Scaling x
        double kz; // Scaling y
        double r;
        double z;
        int size; // matrix size
        double sigma = 100; // S/m
        List<double[,]> Ml; // M locals
        List<double[,]> Gl; // G locals
        double[,] M;
        double[,] G;
        double[,] A; // Matrix
        double[] b; // 'true' decision vector 
        double[] qRe; // finite decision vector (Re)
        double[] qIm; // finite decision vector (Im)
        string q =""; // finite decision vector 
        image[,] portrait; // Matrix portrait
        DrawingGroup drawingGroup = new DrawingGroup();
        List<Element> elements;
        Element field = new Element();
        public MainWindow()
        {
            InitializeComponent();
            buildaxes();
        }

        private void buildaxes()
        {
            GeometryDrawing myGeometryDrawing = new GeometryDrawing();
            GeometryGroup lines = new GeometryGroup();
            myGeometryDrawing.Pen = new Pen(Brushes.Black, 3);
            lines.Children.Add(new LineGeometry(new Point(a, graphImage.Height), new Point(a, 0))); // z
            lines.Children.Add(new LineGeometry(new Point(a, graphImage.Height / 2), new Point(graphImage.Width, graphImage.Height / 2))); // r

            myGeometryDrawing.Geometry = lines;
            drawingGroup.Children.Add(myGeometryDrawing);
            graphImage.Source = new DrawingImage(drawingGroup);

        }

        private void buildMatrix()
        {
            buildlocal();
            buildPortrait();
            // Global M
            size = 2 + Ml.Count * 2;
            M = new double[size, size];
            for (int u = 0; u < Ml.Count; u++)
            {
                for (int i = 0; i < 4; i++)
                {
                    for (int j = 0; j < 4; j++)
                    {
                        M[i + u* 2, j + u * 2] += Ml[u][i, j];
                    }
                }
            }

            // Global G
            G = new double[size, size];
            for (int u = 0; u < Gl.Count; u++)
            {
                for (int i = 0; i < 4; i++)
                {
                    for (int j = 0; j < 4; j++)
                    {
                        G[i + u * 2, j + u * 2] += Gl[u][i, j];
                    }
                }
            }

            // Sum matrixes
            buildA();
            printmat(A, "A");
        }

        // Builds A matrix from M & G
        private void buildA()
        {
            A = new double[size, size];
            for (int i = 0; i < size; i++)
                for (int j = 0;j < size; j++)
                    A[i,j] = -G[i,j] + M[i,j]*(sigma-1);
        }

        // 'True' decision vector
        private void buildAphiRe()
        {
            b = new double[size];
            for (int i = 0; i < size; i++)
                b[i] = 4; // The simplest case, const
        }
        private void buildAphiIm()
        {
            b = new double[size];
            for (int i = 0; i < size; i++)
                b[i] = 4; // The simplest case, const
        }

        private void buildlocal() 
        {
            // p. 233, 5.22 & 5.23
            double[,] Gx = new double[2,2];
            Gx[0, 0] = 1 / elements[0].dr;
            Gx[0, 1] = -1 / elements[0].dr;
            Gx[1, 0] = -1 / elements[0].dr;
            Gx[1, 1] = 1 / elements[0].dr;

            double[,] Gy = new double[2, 2];
            Gy[0, 0] = 1 / elements[0].dz;
            Gy[0, 1] = -1 / elements[0].dz;
            Gy[1, 0] = -1 / elements[0].dz;
            Gy[1, 1] = 1 / elements[0].dz;

            double[,] Mx = new double[2, 2];
            Mx[0, 0] = 2 * elements[0].dr / 6;
            Mx[0, 1] = elements[0].dr/6;
            Mx[1, 0] = elements[0].dr/6;
            Mx[1, 1] = 2 * elements[0].dr/6;

            double[,] My = new double[2, 2];
            My[0, 0] = 2 * elements[0].dz / 6;
            My[0, 1] = elements[0].dz / 6;
            My[1, 0] = elements[0].dz / 6;
            My[1, 1] = 2 * elements[0].dz / 6;

            //G local
            Gl = new List<double[,]>();
            foreach (var element in elements)
            {
                double[,] Gloc = new double[4, 4];
                Gloc[0, 0] = element.lambda * (Gx[0, 0] * My[0,0] + Mx[0, 0] * Gy[0,0]);
                Gloc[0, 1] = element.lambda * (Gx[0, 1] * My[0,0] + Mx[0, 1] * Gy[0,0]);
                Gloc[0, 2] = element.lambda * (Gx[0, 0] * My[0,1] + Mx[0, 0] * Gy[0,1]);
                Gloc[0, 3] = element.lambda * (Gx[0, 1] * My[0, 1] + Mx[0, 1] * Gy[0, 1]);
                Gloc[1, 0] = Gloc[1, 2];
                Gloc[1, 1] = element.lambda * (Gx[1, 1] * My[0, 0] + Mx[1, 1] * Gy[0, 0]);
                Gloc[1, 2] = element.lambda * (Gx[1, 0] * My[0, 1] + Mx[1, 0] * Gy[0, 1]);
                Gloc[1, 3] = element.lambda * (Gx[1, 1] * My[0, 1] + Mx[1, 1] * Gy[0, 1]);
                Gloc[2, 0] = Gloc[0, 2];
                Gloc[2, 1] = Gloc[1, 2];
                Gloc[2, 2] = element.lambda * (Gx[0, 0] * My[1, 1] + Mx[0, 0] * Gy[1, 1]);
                Gloc[2, 3] = element.lambda * (Gx[0, 1] * My[1, 1] + Mx[0, 1] * Gy[1, 1]);
                Gloc[3, 0] = Gloc[0, 3];
                Gloc[3, 1] = Gloc[1, 3];
                Gloc[3, 2] = Gloc[2, 3];
                Gloc[3,3] = element.lambda * (Gx[1, 1] * My[1, 1] + Mx[1, 1] * Gy[1, 1]);
                Gl.Add(Gloc);
            }

            //M local
            Ml = new List <double[,] >();
            foreach (var element in elements)
            {
                double[,] Mloc = new double[4, 4];
                Mloc[0, 0] = element.gamma * Mx[0,0] * My[0,0];
                Mloc[0, 1] = element.gamma * Mx[0,1] * My[0,0];
                Mloc[0, 2] = element.gamma * Mx[0,0] * My[0,1];
                Mloc[0, 3] = element.gamma * Mx[0,1] * My[0,1];
                Mloc[1, 0] = Mloc[0, 1];
                Mloc[1, 1] = element.gamma * Mx[1, 1] * My[0, 0];
                Mloc[1, 2] = element.gamma * Mx[1, 0] * My[0, 1];
                Mloc[1, 3] = element.gamma * Mx[1, 1] * My[0, 1];
                Mloc[2, 0] = Mloc[0, 2];
                Mloc[2, 1] = Mloc[1, 2];
                Mloc[2, 2] = element.gamma * Mx[0, 0] * My[1, 1];
                Mloc[2, 3] = element.gamma * Mx[0, 1] * My[1, 1];
                Mloc[3, 0] = Mloc[0, 3];
                Mloc[3, 1] = Mloc[1, 3];
                Mloc[3, 2] = Mloc[2, 3];
                Mloc[3, 3] = element.gamma * Mx[1, 1] * My[1, 1];
                Ml.Add(Mloc);
            }
        }

        private void Button_Click(object sender, RoutedEventArgs e)
        {
            
            drawingGroup.Children.Clear();
            buildaxes();
            elements = new List<Element>();
            string yval = yVal.Text;
            string xval = xVal.Text;
            if (!(double.TryParse(yval, out double number1)) || !(double.TryParse(xval, out double number2)))
            {
                MessageBox.Show("Ты что вводишь, дуралей?");
            }
            else
            {
                r = double.Parse(xval);
                z = double.Parse(yval);
                // points of the field
                field.p1 = new Point() { X = 0, Y = -z };
                field.p2 = new Point() { X = 0, Y = 0 };
                field.p3 = new Point() { X = r, Y = -z };
                field.p4 = new Point() { X = r, Y = 0 };


                kr = (graphImage.Width - a) / r;
                kz = (graphImage.Height / 2) / z;

                //Draw Field
                GeometryDrawing myGeometryDrawing = new GeometryDrawing();
                GeometryGroup lines = new GeometryGroup();
                myGeometryDrawing.Pen = new Pen(Brushes.Red, 2);

                //Add lines (field borders)
                lines.Children.Add(new LineGeometry(new Point(a, graphImage.Height / 2), new Point(r * kr - a, graphImage.Height / 2))); // up border
                lines.Children.Add(new LineGeometry(new Point(a, graphImage.Height / 2 + z * kz), new Point(r * kr - a, graphImage.Height / 2 -6 + z * kz))); // dowh
                lines.Children.Add(new LineGeometry(new Point(r * kr - a - 6, graphImage.Height / 2), new Point(r * kr - a -6, graphImage.Height / 2 + z * kz))); // right
                lines.Children.Add(new LineGeometry(new Point(a, graphImage.Height / 2), new Point(a, graphImage.Height / 2 + z * kz))); // left


                myGeometryDrawing.Geometry = lines;
                drawingGroup.Children.Add(myGeometryDrawing);
                graphImage.Source = new DrawingImage(drawingGroup);

                MakeElements();

                //Build&Solve
                //Re-part
                buildMatrix();
                buildAphiRe();
                b = solveMatrix(); // Pseudo-decidion
                qRe = new double[size];
                qRe = solveMatrix(); //Real decidion
                //Im-part
                buildAphiIm();
                b = solveMatrix(); // Pseudo-decidion
                qIm = new double[size];
                qIm = solveMatrix(); //Real decidion

                giveDecidion();
            }
        }

        //Completing q-vector
        private void giveDecidion()
        {
            for(int i = 0; i < size; i++)
            {
                q+= qRe[i].ToString() +' ';
                if (qIm[i] != 0)
                {
                    if (qIm[i] > 0) q += "+"; // Positive sign
                    q += qIm[i].ToString() + "i";
                }
                q += '\n';
            }

            //Print
            System.IO.File.WriteAllText("..\\..\\print\\" + "q" + ".txt", q);

        }
        // Print vector or matrix
        private void printvect(double[] vect, string vectname)
        {
            string s = "";
            for (int i = 0; i < vect.Length; i++)
            {
                s += vect[i].ToString() + '\n';
            }

            System.IO.File.WriteAllText("..\\..\\print\\" + vectname  + ".txt", s);
        }
        private void printmat(double[,] mat, string matname)
        {
            int w = Convert.ToInt32(Math.Pow(mat.Length, 0.5));
            string s = "";
            for (int i = 0; i < w; i++)
            {
                for (int j = 0; j < w; j++)
                {
                    s += mat[i, j].ToString() + '\t';
                }
                s += '\n';
            }
            System.IO.File.WriteAllText("..\\..\\print\\" + matname + ".txt", s);
        }

        // Portrait Building
        private void buildPortrait()
        {
            int sizej = int.Parse(xCrush.Text) + 1;
            int sizei = int.Parse(yCrush.Text) + 1;
            portrait = new image[sizei, sizej];
            int num = 0;
            for(int i = 0; i < sizei;i++)
                for (int j = 0;j < sizej; j++)
                {
                    image im = new image();
                    im.n = num;
                    portrait[i, j] = im;
                    num++;
                }
            // Making supporting Matrix with elements numbers. Look at Scheme 1
            int[,] sup = new int[sizei - 1, sizej - 1];
            num = 0; // reuse
            for (int i = 0; i < sizei - 1; i++)
                for (int j = 0; j < sizej - 1; j++)
                {
                    sup[i, j] = num;
                    num++;
                }
            // Search and add elements to corresponding nodes
            for (int i = 0; i < sizei - 1; i++)
                for (int j = 0; j < sizej - 1; j++)
                { 
                    // Define and search nodes
                    int n0 = sup[i, j]  +i;
                    int n1 = n0+1;
                    int n2 = n0+sizej;
                    int n3 = n2 + 1;
                    // Go all nodes 
                    for (int k = 0; k < sizei; k++)
                        for (int m = 0; m < sizej; m++)
                        {
                            if (portrait[k,m].n == n0 || portrait[k, m].n == n1 || 
                                portrait[k, m].n == n2 || portrait[k, m].n == n3)
                                portrait[k, m].els.Add(sup[i,j]);
                        }

                }
                }

        // building of elements
        private void MakeElements()
        {
            string zcrush = yCrush.Text;
            string xcrush = xCrush.Text;
            if (!(double.TryParse(zcrush, out double number1)) || !(double.TryParse(xcrush, out double number2)))
            {
                MessageBox.Show("Ты что вводишь, дуралей?");
            }
            else if (double.Parse(xcrush) <= 0 || double.Parse(zcrush) <= 0)
            {
                MessageBox.Show("Ты что вводишь, дуралей?");
            }
            else
            {
                double dr = r / double.Parse(xcrush);
                double dz = z / double.Parse(zcrush);

                for (int i = 1; i <= double.Parse(xcrush); i++)
                    for (int j = 1; j <= double.Parse(zcrush); j++)
                    {
                        Element el = new Element();
                        el.p1 = new Point() { X = (i - 1) * dr, Y = -j * dz };
                        el.p2 = new Point() { X = (i - 1) * dr, Y = -(j - 1) * dz };
                        el.p3 = new Point() { X = i * dr, Y = -j * dz };
                        el.p4 = new Point() { X = i * dr, Y = -(j - 1) * dz };

                        el.dr = dr;
                        el.dz = dz;

                        elements.Add(el);
                    }

                // Paint
                GeometryDrawing myGeometryDrawing = new GeometryDrawing();
                GeometryGroup lines = new GeometryGroup();
                myGeometryDrawing.Pen = new Pen(Brushes.Red, 3);
                for (int i = 1; i < int.Parse(xcrush); i++)
                {
                    lines.Children.Add(new LineGeometry(new Point(i * dr * kr, (graphImage.Height -6)/ 2), new Point(i * dr * kr, (graphImage.Height - 6) / 2 + z * kz)));
                }
                for (int j = 1; j < int.Parse(zcrush); j++)
                {
                    lines.Children.Add(new LineGeometry(new Point(a, j * dz * kz + (graphImage.Height -6)/ 2), new Point(r * kr - a, j * dz * kz + (graphImage.Height - 6) / 2)));
                }
                myGeometryDrawing.Geometry = lines;
                drawingGroup.Children.Add(myGeometryDrawing);
                graphImage.Source = new DrawingImage(drawingGroup);
            }
        }


        // Gauss solver
        private double[] solveMatrix()
        {
            double[] ans = new double[b.Length];
            solveSLAE(ans);
            return ans;
        }
        // Solves SOLaE
        public void solveSLAE(double[] ans)
        {
            int nSLAE = b.Length;
            if (ans.Length != nSLAE)
                throw new Exception("Size of the input array is not compatable with size of SLAE");




            for (int i = 0; i < nSLAE; i++)
            {
                double del = A[i, i];
                double absDel = Math.Abs(del);
                int iSwap = i;


                for (int j = i + 1; j < nSLAE; j++) // ищем максимальный элемент по столбцу
                {
                    if (absDel < Math.Abs(A[j, i]))
                    {
                        del = A[j, i];
                        absDel = Math.Abs(del);
                        iSwap = j;
                    }
                }

                if (iSwap != i)
                {
                    double buf;
                    for (int j = i; j < nSLAE; j++)
                    {
                        buf = A[i, j];
                        A[i, j] = A[iSwap, j];
                        A[iSwap, j] = buf;
                    }
                    buf = b[i];
                    b[i] = b[iSwap];
                    b[iSwap] = buf;
                }

                for (int j = i; j < nSLAE; j++)
                    A[i, j] /= del;

                b[i] /= del;

                for (int j = i + 1; j < nSLAE; j++)
                {
                    if (A[j, i] == 0) continue;

                    double el = A[j, i];
                    for (int k = i; k < nSLAE; k++)
                    {
                        A[j, k] -= A[i, k] * el;
                    }

                    b[j] -= b[i] * el;
                }
            }

            for (int i = nSLAE - 1; i > -1; i--)
            {
                for (int j = i + 1; j < nSLAE; j++)
                    b[i] -= ans[j] * A[i, j];
                ans[i] = b[i];
            }
        }


        

    }


    // The struct of an element
    /*
     p2-----------p4
     |             |
     |             |
     |             |
     |             |
     p1-----------p3
         
     Height: dz
     Width: dr
     
     */
    public class Element
    {
        public Point p1 { get; set; }
        public Point p2 { get; set; }
        public Point p3 { get; set; }
        public Point p4 { get; set; }

        public double dr;

        public double dz;

        public double lambda = 1 / (1e-7 * 4 * Math.PI);// 1/myu0
        public double gamma = 1 ;// ??
    }

    public class image
    {
        public int n;
        public List<int> els = new List<int>();
    }
}
