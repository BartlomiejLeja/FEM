using System;
using FEM.FemGrid;

namespace FEM
{
    class Program
    {
        static void Main(string[] args)
        {
           var cofficentsProvider = new CofficientsProvider();
           var integral = new Integral();
           var result2p = integral.Calculate(cofficentsProvider.coefficients2p);
           var result3p = integral.Calculate(cofficentsProvider.coefficients3p);
           Console.WriteLine("Double integrating");
           Console.WriteLine($"Result of two point double integrating is {result2p}");
           Console.WriteLine($"Result of three point double integrating is {result3p}");
           Console.WriteLine($"Difference is {result3p - result2p}");

           Console.WriteLine("Fem grid");
           var grid = new Grid();
           grid.FillFemGrid();
            grid.CalculateEveryThing();

           // var universalElement = new UniversalElement();
            //double[,] pointTab = new Double[4, 2]
            //{
            //    {0,0 },
            //    {2,0 },
            //    {2,2 },
            //    {0,2 }
            //};

            //universalElement.CalculateJackobian(pointTab, 1);

            Console.ReadKey();
        }
    }
}
