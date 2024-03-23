using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography;
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
using static System.Net.Mime.MediaTypeNames;

namespace Electric
{
    /// <summary>
    /// Логика взаимодействия для MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        double sr; // r for source
        double sz; // z for source
        double a = 3; // Only for drawning
        double kr; // Scaling x
        double kz; // Scaling y
        double r;
        double z;
        double w = 1;
        double dr; // element's width
        double dz; // element's deepth
        int sizei; // matrix size
        int sizej; // matrix size
        int size; // Global matrix size
        bool iflag = false;
        double sigma = 100; // S/m
        List<double[,]> Ml; // M locals
        List<double[,]> Gl; // G locals
        double[,] M;
        double[,] G;
        double[,] A; // Matrix
        double[,] Ac; // cos matrix (Re)
        double[,] As; // sin Matrix (Im)
        double[] b; // 'true' decision vector 
        double[] q; // finite decision vector (Re)
        int[,] Nods;
        //string q =""; // finite decision vector 
        img[,] portrait; // Matrix portrait
        DrawingGroup drawingGroup = new DrawingGroup();
        List<Element> elements;
        Element field = new Element();
        Source source = new Source();
        
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
            size = (sizei + 1) * (sizej + 1);
            M = new double[size, size]; // Size = number of nodes * number of nodes
            for (int i = 0; i < sizej; i++)
                for (int j = 0; j < sizei; j++)
                {
                    // Go local Matrix
                    int num = portrait[i, j].n;
                    for(int k = 0; k < 4; k++)
                        for(int l = 0; l < 4; l++)
                        {
                            // Global indexes
                            int ig = portrait[i, j].els[k];
                            int jg = portrait[i, j].els[l];
                            M[ig, jg] += Ml[num][k,l];
                        }
                }

            // Global G
            G = new double[size, size];
            for (int i = 0; i < sizej; i++)
                for (int j = 0; j < sizei; j++)
                {
                    // Go local Matrix
                    int num = portrait[i, j].n;
                    for (int k = 0; k < 4; k++)
                        for (int l = 0; l < 4; l++)
                        {
                            // Global indexes
                            int ig = portrait[i, j].els[k];
                            int jg = portrait[i, j].els[l];
                            G[ig, jg] += Gl[num][k, l];
                        }
                }

            // Sum matrixes
            buildA();
            printmat(A, "A");
            printmat(M, "M");
            printmat(G, "G");
        }

        // Builds A matrix from M & G

        private void buildsincos() //As & Ac
        {
            A = new double[size * 2, size * 2];
            Ac = new double[size, size];
            As = new double[size, size];
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                {
                    /*Ac[i, j] = -G[i, j] + M[i, j] * (sigma - 1);
                    As[i, j] = -G[i, j] - M[i, j] * (sigma + 1);*/
                    Ac[i, j] = -G[i, j] + M[i, j] * (sigma - w)*w;
                    As[i, j] = -G[i, j] - M[i, j] * (w*sigma + w*w);
                }
        }
        private void buildA()
        {
            buildsincos();
            // Look at Scheme 2
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                {
                    int i0 = i * 2; 
                    int j0= j*2;
                    A[i0, j0]= As[i,j];
                    A[i0, j0 + 1]= -Ac[i,j];
                    A[i0 + 1, j0]= Ac[i,j];
                    A[i0 + 1, j0 + 1]= As[i,j];
                }



        }

        // 'True' decision vector
        private void buildAphi()
        {
            for (int i = 0; i < size*2; i++)
                if(i%2==0)
                b[i] = 4; // The simplest case, const
            else b[i] = 3;
        }

        //testing
        private void buildAtest()
        {
            int n = 0; // Number of a node
            
            buildsincos();
            
            buildA();
            // Inserting requred values
            for (int i = 0; i <= sizej; i++)
                for (int j = 0; j <= sizei; j++)
                {
                    if (i == 0 || j == 0 || i == sizei || j == sizej)
                    {
                        n = Nods[i, j]*2;
                        // Changing strings in A
                        for(int k = 0; k < size*2; k++)
                        {
                            A[n,k] = 0;
                            A[n+1,k] = 0;
                        }
                        A[n, n] = 1;
                        A[n + 1, n + 1] = 1;
                    }

                }
            
        }

        // Changes b-vector
        private void testbconst() 
        {
            for (int i = 0; i <= sizej; i++)
                for (int j = 0; j <= sizei; j++)
                {
                    if (i == 0 || j == 0 || i == sizei || j == sizej)
                    {
                        int n = Nods[i, j];
                        b[n * 2] = 7;
                        b[n * 2 + 1] = 7;
                    }
                }
        }
        private void testbr(int pow) 
        {
            double drr = double.Parse(xVal.Text) / sizei;
            double drz = double.Parse(yVal.Text) / sizej;
            for (int i = 0; i <= sizej; i++)
                for (int j = 0; j <= sizei; j++)
                {
                    if (i == 0 || j == 0 || i == sizei || j == sizej)
                    {
                        double r = Math.Sqrt(Math.Pow(drr*j, 2) + Math.Pow(drz*i, 2));
                        int n = Nods[i, j];
                        b[n * 2] = Math.Pow(r, pow);
                        b[n * 2 + 1] = Math.Pow(r, pow);
                    }
                }
        }
        private void test()
        {
            double[] qtest = new double[size];

            // Get right part
            buildAphi();
            buildA();
            b = solveMatrix();
            printvect(b, "btest");

            // Get A for test
            buildAtest();
            printmat(A, "Atest");

            // Change right part
            testbconst(); 
            printvect(b, "btestconst");

            // Get result
            qtest = solveMatrix();
            printvect(qtest, "qtest");

            // f = r
            buildAphi();
            buildA();
            b = solveMatrix();
            // Get A for test
            buildAtest();
            // Change right part
            testbr(1);
            printvect(b, "btestr");
            qtest = solveMatrix();
            printvect(qtest, "qtestr");

            // f = r^2
            buildAphi();
            buildA();
            b = solveMatrix();
            // Get A for test
            buildAtest();
            // Change right part
            testbr(2);
            printvect(b, "btestr^2");
            qtest = solveMatrix();
            printvect(qtest, "qtestr^2");

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
                Gloc[1, 0] = Gloc[0, 1];
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
                // All coefficient for M
                double koef = element.sigma - element.lambda * 1/element.r;
                double[,] Mloc = new double[4, 4];
                Mloc[0, 0] = element.gamma * Mx[0,0] * My[0,0] * koef;
                Mloc[0, 1] = element.gamma * Mx[0,1] * My[0,0] * koef;
                Mloc[0, 2] = element.gamma * Mx[0,0] * My[0,1] * koef;
                Mloc[0, 3] = element.gamma * Mx[0, 1] * My[0, 1] * koef;
                Mloc[1, 0] = Mloc[0, 1];
                Mloc[1, 1] = element.gamma * Mx[1, 1] * My[0, 0] * koef;
                Mloc[1, 2] = element.gamma * Mx[1, 0] * My[0, 1] * koef;
                Mloc[1, 3] = element.gamma * Mx[1, 1] * My[0, 1] * koef;
                Mloc[2, 0] = Mloc[0, 2];
                Mloc[2, 1] = Mloc[1, 2];
                Mloc[2, 2] = element.gamma * Mx[0, 0] * My[1, 1] * koef;
                Mloc[2, 3] = element.gamma * Mx[0, 1] * My[1, 1] * koef;
                Mloc[3, 0] = Mloc[0, 3];
                Mloc[3, 1] = Mloc[1, 3];
                Mloc[3, 2] = Mloc[2, 3];
                Mloc[3, 3] = element.gamma * Mx[1, 1] * My[1, 1] * koef;
                Ml.Add(Mloc);
            }
        }

        private void Button_Click(object sender, RoutedEventArgs e)
        {
            //Source
            if (!(double.TryParse(Sourcer.Text, out double num1)) || !(double.TryParse(Sourcez.Text, out double num2)))
            {
                MessageBox.Show("Wrong!");
                return;
            }
            else
            {
             sr = double.Parse(Sourcer.Text);
             sz = double.Parse(Sourcez.Text);
                if(sz < 0 || sr < 0)
                {
                    MessageBox.Show("Wrong!");
                    return;
                }
            }
            source.rz = new Point(sr, -sz);
            // Elements
            sigma = double.Parse(sigVal.Text);
            w = double.Parse(wVal.Text);
            drawingGroup.Children.Clear();
            buildaxes();
            elements = new List<Element>();
            string yval = yVal.Text;
            string xval = xVal.Text;
            if (!(double.TryParse(yval, out double number1)) || !(double.TryParse(xval, out double number2)))
            {
                MessageBox.Show("Wrong!");
                return;
            }
                r = double.Parse(xval);
                z = double.Parse(yval);
                // points of the field
                field.p1 = new Point() { X = 0, Y = -z };
                field.p2 = new Point() { X = 0, Y = 0 };
                field.p3 = new Point() { X = r, Y = -z };
                field.p4 = new Point() { X = r, Y = 0 };


                kr = (graphImage.Width - a) / r;
                kz = ((graphImage.Height - a) / 2) / z;

                //Draw Field
                GeometryDrawing myGeometryDrawing = new GeometryDrawing();
                GeometryGroup lines = new GeometryGroup();
                myGeometryDrawing.Pen = new Pen(Brushes.Red, 1);

                //Add lines (field borders)
                lines.Children.Add(new LineGeometry(new Point(a, graphImage.Height / 2), new Point(r * kr - a, graphImage.Height / 2))); // up border
                lines.Children.Add(new LineGeometry(new Point(a, graphImage.Height / 2 + z * kz), new Point(r * kr - a, graphImage.Height / 2 + z * kz))); // down
                lines.Children.Add(new LineGeometry(new Point(r * kr - a, graphImage.Height / 2), new Point(r * kr - a, graphImage.Height / 2 + z * kz))); // right
                lines.Children.Add(new LineGeometry(new Point(a, graphImage.Height / 2), new Point(a, graphImage.Height / 2 + z * kz))); // left


                myGeometryDrawing.Geometry = lines;
                drawingGroup.Children.Add(myGeometryDrawing);
                graphImage.Source = new DrawingImage(drawingGroup);

                MakeElements();

                 // Find source's node number
                SetSourcen();
                double dphi = CountPotentials();
                //Build&Solve
                //Re-part
                buildMatrix();
                b = new double[size * 2];
              //  b[source.n] = 1;
                buildAphi();
                b = solveMatrix(); // Pseudo-decidion
                printvect(b, "b");
                q = new double[size*2];
                q = solveMatrix(); //Real decidion
                giveDecidion();
                test();
            
        }

        // Func for source node number
        private void SetSourcen()
        {
            int currn = 0;
            int redge = Convert.ToInt32(xCrush.Text)*2;
            int drparts = redge + 1;
            //int zparts = Convert.ToInt32(yCrush.Text)*2;
            foreach(var el in elements)
            {
                double diff = Math.Abs(Math.Abs(el.p2.X) - Math.Abs(source.rz.X));
                diff += Math.Abs(Math.Abs(el.p2.Y) - Math.Abs(source.rz.Y));
                if (diff < 0.0001)
                {
                    source.n = currn;
                    break;
                }
                currn++;
                // for right edge elements
                if (currn % redge == 0)
                {
                    currn++;
                    redge += drparts;
                }
            }
        }

        //Completing q-vector
        private void giveDecidion()
        {
            string qv = "";  // q-vector
           for(int i = 0; i < q.Length; i++)
            {
                qv+= q[i].ToString() +' ';
                if (iflag)
                {
                    if (i % 2 == 1 && q[i] != 0)
                    {
                        qv += 'i';
                    }
                }
                qv += '\n';
            }

            //Print
            System.IO.File.WriteAllText("..\\..\\print\\" + "q" + ".txt", qv);

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
            int l = Convert.ToInt32(Math.Pow(mat.Length, 0.5));
            string s = "";
            for (int i = 0; i < l; i++)
            {
                for (int j = 0; j < l; j++)
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
            
            portrait = new img[sizej, sizei];
            int num = 0;
            for(int i = 0; i < sizej;i++)
                for (int j = 0;j < sizei; j++)
                {
                    img im = new img();
                    im.n = num;
                    portrait[i, j] = im;
                    num++;
                }
            
            // Search and add elements to corresponding nodes
            for (int i = 0; i < sizej; i++)
                for (int j = 0; j < sizei; j++)
                { 
                    // Define nodes
                    int n2 = portrait[i, j].n + i; // 2
                    int n3 = n2+1; // 3
                    int n0 = n2+(sizei +1); // 0
                    int n1 = n0 + 1; // 1

                    // Add nodes in the order
                    portrait[i, j].els.Add(n0);
                    portrait[i, j].els.Add(n1);
                    portrait[i, j].els.Add(n2);
                    portrait[i, j].els.Add(n3);
                }
                }

        // building of elements
        private void MakeElements()
        {
            double rdis = Convert.ToDouble(rdischarge.Text);
            double zdis = Convert.ToDouble(zdischarge.Text);
            int insr = Convert.ToInt32(xCrush.Text);                              // initial step r
            int insz = Convert.ToInt32(yCrush.Text);                              // initial step z
            Point zeroPoint = new Point(Convert.ToDouble(Sourcer.Text) * kr +a,   // for drawing
                graphImage.Height / 2 + Convert.ToDouble(Sourcez.Text) * kz);
            Point source = new Point(Convert.ToDouble(Sourcer.Text),              // for counting
                -Convert.ToDouble(Sourcez.Text));
            // Maing supporting areas, look Scheme3
            List<Element> supp1 = new List<Element>();
            List<Element> supp2 = new List<Element>();
            List<Element> supp3 = new List<Element>();
            List<Element> supp4 = new List<Element>();

            // Making elements
            // Counting first drl &
            double ksumx = 0;
            for (int i = 0; i < Convert.ToInt32(xCrush.Text); i++)
                ksumx += Math.Pow(rdis, i);
            double drl = source.X / ksumx;
            double ksumy = 0;
            for (int i = 0; i < Convert.ToInt32(yCrush.Text); i++)
                ksumy += Math.Pow(zdis, i);
            //
            double dzu = source.Y / ksumy;                                   // from up, < 0
            double dzd = (-Convert.ToDouble(yVal.Text) - source.Y)/ ksumy;   // from down, < 0
            //double drl = source.X / Convert.ToDouble(xCrush.Text);                                   // from left
            double drr =(Convert.ToDouble(xVal.Text) - source.X) / ksumx;    // from right
            //prevoius x & y lines, look Scheme 3
            double ybottom = -Convert.ToDouble(yVal.Text); 
            double xrightside = Convert.ToDouble(xVal.Text);

            Point p3 = source;
            Point p4 = new Point() { X = source.X, Y = source.Y + insz};
            Point p2 = new Point() { X = source.X - insr, Y = source.Y +insz};
            Point p1 = new Point() { X = source.X - insr, Y = source.Y };
            double rst = insr;
            double zst = insz;
            int rpartsleft;
            int zpartsup = 0;
            double rd, zd;
            // Making elements in the first area
            while (p2.Y <= 0)
            {  
                while (p1.X >= 0 )
                {
                    Element el1 = new Element();
                    el1.p1 = p1;
                    el1.p2 = p2;
                    el1.p4 = p4;
                    el1.p3 = p3;

                    el1.dr = p4.X - p2.X;
                    el1.dz = p4.Y - p3.Y;
                    rd = Math.Abs(source.X - p2.X) - 0.5*el1.dr;
                    zd = Math.Abs(source.Y -  p2.Y) - 0.5*el1.dz; 
                    el1.r = Math.Sqrt(Math.Pow(rd, 2) + Math.Pow(zd, 2));

                    supp1.Add(el1);
                    // itteration
                    rst *= rdis;
                    p3 = p1;
                    p4 = p2;
                    p1.X -= rst;
                    p2.X -= rst;
                }
                // left border
                int last = supp1.Count - 1;
                if (supp1[last].p1.X > 0)
                {
                    var pb = supp1[last].p1;
                    pb.X = 0;
                    supp1[last].p1 = pb;
                    pb = supp1[last].p2;
                    pb.X = 0;
                    supp1[last].p2= pb;
                }
                // next row
                rst = insr;
                p3.X = source.X;
                p3.Y += zst;
                p1 = p3;
                p1.X -= rst;
                zst *= zdis;
                p4 = p3;
                p4.Y += zst;
                p2 = p4;
                p2.X -= rst;
                zpartsup++;
            }
            // up border
            rpartsleft = supp1.Count / zpartsup;
            if (supp1[supp1.Count - 1].p2.Y < 0)
            {
                for (int i = 0; i < rpartsleft; i++)
                {
                    int ind = supp1.Count - 1 - i;
                    var pb = supp1[ind].p2;
                    pb.Y = 0;
                    supp1[ind].p2 = pb;
                    pb = supp1[ind].p4;
                    pb.Y = 0;
                    supp1[ind].p4 = pb;
                }
            }

            // Making elements in the second area
            p4 = source;
            p3 = new Point() { X = source.X, Y = source.Y - insz };
            p1 = new Point() { X = source.X - insr, Y = source.Y - insz };
            p2 = new Point() { X = source.X - insr, Y = source.Y};
            zst = insz;
            rst = insr;
            int zpartsdown = 0;

            while (p1.Y >= -Convert.ToDouble(yVal.Text))
            {
                while (p1.X >= 0)
                {
                    Element el1 = new Element();
                    el1.p1 = p1;
                    el1.p2 = p2;
                    el1.p4 = p4;
                    el1.p3 = p3;

                    el1.dr = p4.X - p2.X;
                    el1.dz = p4.Y - p3.Y;
                    rd = Math.Abs(source.X - p2.X) - 0.5 * el1.dr;
                    zd = Math.Abs(source.Y - p2.Y) + 0.5 * el1.dz;
                    el1.r = Math.Sqrt(Math.Pow(rd, 2) + Math.Pow(zd, 2));

                    supp2.Add(el1);
                    // itteration
                    rst *= rdis;
                    p3 = p1;
                    p4 = p2;
                    p1.X -= rst;
                    p2.X -= rst;
                    el1.dr = 1;
                }
                // left border
                int last = supp2.Count - 1;
                if (supp2[last].p1.X > 0)
                {
                    var pb = supp2[last].p1;
                    pb.X = 0;
                    supp2[last].p1 = pb;
                    pb = supp2[last].p2;
                    pb.X = 0;
                    supp2[last].p2 = pb;
                }
                // next row
                rst = insr;
                p4.X = source.X;
                p4.Y -= zst;
                p2 = p4;
                p2.X -= rst;
                zst *= zdis;
                p3 = p4;
                p3.Y -= zst;
                p1 = p3;
                p1.X -= rst;
                zpartsdown++;
            }
            // down border
            if (supp2[supp2.Count - 1].p1.Y > -Convert.ToDouble(yVal.Text))
            {
                for (int i = 0; i < rpartsleft; i++)
                {
                    int ind = supp2.Count - 1 - i;
                    var pb = supp2[ind].p1;
                    pb.Y = -Convert.ToDouble(yVal.Text);
                    supp2[ind].p1 = pb;
                    pb = supp2[ind].p3;
                    pb.Y = -Convert.ToDouble(yVal.Text);
                    supp2[ind].p3 = pb;
                }
            }

            // Making elements in the third area
            p1 = source;
            p4 = new Point() { X = source.X + insr, Y = source.Y + insz };
            p2 = new Point() { X = source.X , Y = source.Y + insz };
            p3 = new Point() { X = source.X + insr, Y = source.Y };
            zst = insz;
            rst = insr;
            int rpartsright;

            while (p2.Y <= 0)
            {
                while (p4.X <= Convert.ToDouble(xVal.Text))
                {
                    Element el1 = new Element();
                    el1.p1 = p1;
                    el1.p2 = p2;
                    el1.p4 = p4;
                    el1.p3 = p3;

                    el1.dr = p4.X - p2.X;
                    el1.dz = p4.Y - p3.Y;
                    rd = Math.Abs(source.X - p2.X) + 0.5 * el1.dr;
                    zd = Math.Abs(source.Y - p2.Y) - 0.5 * el1.dz;
                    el1.r = Math.Sqrt(Math.Pow(rd, 2) + Math.Pow(zd, 2));

                    supp3.Add(el1);
                    // itteration
                    rst *= rdis;
                    p1 = p3;
                    p2 = p4;
                    p3.X += rst;
                    p4.X += rst;
                    el1.dr = 1;
                }
                // right border
                int last = supp3.Count - 1;
                if (supp3[last].p4.X < Convert.ToDouble(xVal.Text))
                {
                    var pb = supp3[last].p4;
                    pb.X = Convert.ToDouble(xVal.Text);
                    supp3[last].p4 = pb;
                    pb = supp3[last].p3;
                    pb.X = Convert.ToDouble(xVal.Text);
                    supp3[last].p3 = pb;
                }
                // next row
                rst = insr;
                p1.X = source.X;
                p1.Y += zst;
                p3 = p1;
                p3.X += rst;
                zst *= zdis;
                p2 = p1;
                p2.Y += zst;
                p4 = p2;
                p4.X += rst;
            }
            rpartsright = supp3.Count / zpartsup;
            // up border
            if (supp3[supp3.Count - 1].p2.Y < 0)
            {
                for (int i = 0; i < rpartsleft; i++)
                {
                    int ind = supp3.Count - 1 - i;
                    var pb = supp3[ind].p2;
                    pb.Y = 0;
                    supp3[ind].p2 = pb;
                    pb = supp3[ind].p4;
                    pb.Y = 0;
                    supp3[ind].p4 = pb;
                }
            }

            // Making elements in the fourth area
            p2 = source;
            p4 = new Point() { X = source.X + insr, Y = source.Y};
            p1 = new Point() { X = source.X, Y = source.Y - insz };
            p3 = new Point() { X = source.X + insr, Y = source.Y - insz };
            zst = insz;
            rst = insr;

            while (p1.Y >= -Convert.ToDouble(yVal.Text))
            {
                while (p3.X <= Convert.ToDouble(xVal.Text))
                {
                    Element el1 = new Element();
                    el1.p1 = p1;
                    el1.p2 = p2;
                    el1.p4 = p4;
                    el1.p3 = p3;

                    el1.dr = p4.X - p2.X;
                    el1.dz = p4.Y - p3.Y;
                    rd = Math.Abs(source.X - p2.X) + 0.5 * el1.dr;
                    zd = Math.Abs(source.Y - p2.Y) + 0.5 * el1.dz;
                    el1.r = Math.Sqrt(Math.Pow(rd, 2) + Math.Pow(zd, 2));

                    supp4.Add(el1);
                    // itteration
                    rst *= rdis;
                    p1 = p3;
                    p2 = p4;
                    p3.X += rst;
                    p4.X += rst;
                    el1.dr = 1;
                }
                // right border
                int last = supp4.Count - 1;
                if (supp4[last].p3.X < Convert.ToDouble(xVal.Text))
                {
                    var pb = supp4[last].p3;
                    pb.X = Convert.ToDouble(xVal.Text);
                    supp4[last].p3 = pb;
                    pb = supp4[last].p4;
                    pb.X = Convert.ToDouble(xVal.Text);
                    supp4[last].p4 = pb;
                }
                // next row
                rst = insr;
                p2.X = source.X;
                p2.Y -= zst;
                zst *= zdis;
                p4 = p2;
                p4.X += rst;
                p1 = p2;
                p1.Y -= zst;
                p3 = p4;
                p3.Y -= zst;
            }
            // down border
            if (supp4[supp4.Count - 1].p1.Y > -Convert.ToDouble(yVal.Text))
            {
                for (int i = 0; i < rpartsleft; i++)
                {
                    int ind = supp4.Count - 1 - i;
                    var pb = supp4[ind].p1;
                    pb.Y = -Convert.ToDouble(yVal.Text);
                    supp4[ind].p1 = pb;
                    pb = supp4[ind].p3;
                    pb.Y = -Convert.ToDouble(yVal.Text);
                    supp4[ind].p3 = pb;
                }
            }

            // Seting elements in the right order
            for (int j = zpartsup - 1; j >=0  ; j--)
            {
                for (int i = rpartsleft - 1; i >= 0; i--)
                {
                    elements.Add(supp1[i + j* rpartsleft]);
                }
                for (int i = 0; i < rpartsright; i++)
                {
                    elements.Add(supp3[i + j * rpartsright]);
                }
            }
            for (int j = 0; j < zpartsdown ; j++)
            {
                for (int i = rpartsleft - 1; i >= 0; i--)
                {
                    elements.Add(supp2[i + j * rpartsleft]);
                }
                for (int i = 0; i < rpartsright; i++)
                {
                    elements.Add(supp4[i + j * rpartsright]);
                }
            }

            //sizes
            // for matrixes
            sizei = rpartsleft + rpartsright;
            sizej = zpartsdown+zpartsup;

            // Paint
            GeometryDrawing myGeometryDrawing = new GeometryDrawing();
                    GeometryDrawing sourceDrawing = new GeometryDrawing();
                GeometryGroup lines = new GeometryGroup();
                myGeometryDrawing.Pen = new Pen(Brushes.Red, 1);
                sourceDrawing.Pen = new Pen(Brushes.Blue, 1);
            foreach(var e in elements)
            {
                lines.Children.Add(new LineGeometry(new Point(e.p1.X*kr+a, -e.p1.Y*kz + graphImage.Height/2),
                    new Point(e.p2.X * kr + a, -e.p2.Y * kz + graphImage.Height / 2)));

                lines.Children.Add(new LineGeometry(new Point(e.p4.X * kr + a, -e.p4.Y * kz + graphImage.Height / 2),
                    new Point(e.p3.X * kr +a, -e.p3.Y * kz + graphImage.Height / 2)));

                lines.Children.Add(new LineGeometry(new Point(e.p2.X * kr + a, -e.p2.Y * kz + graphImage.Height / 2),
                    new Point(e.p4.X * kr + a, -e.p4.Y * kz + graphImage.Height / 2)));

                lines.Children.Add(new LineGeometry(new Point(e.p1.X * kr +a, -e.p1.Y * kz + graphImage.Height / 2),
                    new Point(e.p3.X * kr +a, -e.p3.Y * kz + graphImage.Height / 2)));
            }
                drawingGroup.Children.Add(myGeometryDrawing);
                // Paint source
                EllipseGeometry ellipse = new EllipseGeometry();
            ellipse.Center = zeroPoint;
            ellipse.RadiusX = 3;
            ellipse.RadiusY = 3;

            myGeometryDrawing.Geometry = lines;
            sourceDrawing.Geometry = ellipse;
            drawingGroup.Children.Add(sourceDrawing);
            graphImage.Source = new DrawingImage(drawingGroup);
            
            // Matrixes sizes
            /*sizej = int.Parse(xCrush.Text) * 2;
            sizei = int.Parse(yCrush.Text) * 2;*/
            int n = 0;
            Nods = new int[sizej + 1, sizei + 1];
            for (int i = 0; i <= sizej; i++)
                for (int j = 0; j <= sizei; j++)
                {
                    Nods[i, j] = n;
                    n++;
                }
        }

        // Forming source
        private double CountPotentials()
        {
            // Physical values
            double J = 0.7;                                         // A
            double f = 1;                                           // MHz
            double myu = 4 * Math.PI * Math.Pow(10, -7);            // H/m
            double L1 = -source.rz.Y;                               // m
            double L2 = source.rz.Y + Convert.ToDouble(yVal.Text);  // m
            double ro = 1000;                                       // kg/m^3

            // Coils
            Coil I1 = new Coil() { n = 100, S = 0.002};
            Coil I2 = new Coil() { n = 200, S = 0.003 };
            Coil G = new Coil() { n = 150, S = 0.001 };
            I1.M = I1.S*I1.n;
            I2.M = I2.S*I2.n;
            G.M = G.S*G.n*J;

            // Count
            double p1 = Math.Sqrt(Math.PI * f * myu / ro) * L1;
            double dL = L1 - L2;
            double deltaL = dL / L1;
            double dphi = p1 * deltaL;
            dphi -= Math.Atan(p1 * deltaL / (1 + p1*(2 - deltaL) + 2*p1*p1*(1 - deltaL)));
            return dphi;
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

        private void graphImage_MouseMove(object sender, MouseEventArgs e)
        {
            double mr = e.GetPosition(graphImage).X;
            double mz = e.GetPosition(graphImage).Y;
            double nr;
            double nz;
            for (int i = 0; i <= sizei; i++)
                for (int j = 0; j <= sizej; j++)
                
                {
                     nr = dr * kr * j;      ///kx^...
                     nz = graphImage.Height / 2 + dz * kz * i;
                    if (nr - 3 < mr && nr + 3 > mr)
                    {
                        if (nz - 3 < mz && nz + 3 > mz)
                        {
                            ToolTip tooltip = new ToolTip { Content = Nods[i, j].ToString() };
                            graphImage.ToolTip = tooltip;
                        }
                    }
                    
                }

        }

        private void CheckBox_Checked(object sender, RoutedEventArgs e)
        {
            iflag = !iflag;
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
        public double sigma = 1;
        public double r;
    }

    public class img
    {
        public int n;
        public List<int> els = new List<int>();
        
    }
    public class Source
    {
        public Point rz;
        public int n;
    }
    public class Coil
    {
        public int n;
        public double S;
        public double M;
    }
}
