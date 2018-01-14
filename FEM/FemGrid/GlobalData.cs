using System;
using System.Globalization;
using System.IO;
using System.Text;

namespace FEM.FemGrid
{
    class GlobalData
    {
        public static double H { get; set; }
        public static double  B { get; set; }
        public static int nH { get; set; }
        public static int nB { get; set; }
        public static int nn { get; set; }
        public static int ne { get; set; }
        public static int k { get; set; }//wsp. ciepła
        public static int alfa { get; set; }
        public static int c { get; set; }
        public static int ro { get; set; }
        public static int temp { get; set; }
        public static double L1 { get; set; }
        public static double L2 { get; set; }

        public static void GetGlobalData(string path)
        {
            string[] lines = File.ReadAllLines(path, Encoding.UTF8);

            H = 0.1; /*double.Parse(lines[0].Split(' ')[0], CultureInfo.InvariantCulture);*/
            B = 0.1; /*double.Parse(lines[1].Split(' ')[0], CultureInfo.InvariantCulture);*/
            nH = 4; /*Convert.ToInt16( lines[2].Split(' ')[0]);*/
            nB = 4; /*Convert.ToInt16(lines[3].Split(' ')[0]);*/
            k = 25; /*Convert.ToInt16(lines[4].Split(' ')[0]);*/
            alfa = 300;/* Convert.ToInt16(lines[5].Split(' ')[0]);*/
            c = 700; /*Convert.ToInt16(lines[6].Split(' ')[0]);*/
            ro = 7800; /*Convert.ToInt16(lines[7].Split(' ')[0]);*/
            temp = 1200; /*Convert.ToInt16(lines[7].Split(' ')[0]);*/

            nn = nH * nB;
            ne = (nH-1) * (nB-1);
            L1 = B / (nB - 1);
            L2 = H / (nH - 1);
        }
    }
}
