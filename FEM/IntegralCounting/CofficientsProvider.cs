namespace FEM
{
    class CofficientsProvider
    {
        public double[,] coefficients2p = new double[,]
           {
                {-0.577, 1 },
                {0.577,1 }
           };

        public double[,] coefficients3p = new double[,]
         {
                 {-0.7745, 5.0/9.0 },
                 {0,8.0/9.0 },
                 {0.7745,5.0/9.0 }
         };
    }
}
