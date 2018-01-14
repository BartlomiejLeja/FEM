namespace FEM
{
    class Integral
    {
        private double function(double x, double y)
        {
            return 2 * x*x * y*y + 6*x + 5;
        }
 
        public double Calculate(double[,] coefficients)
        {
            double result=0;
            int npc = coefficients.Length / 2;

            for (int i=0;i<npc;i++)
            {
                for (int j = 0; j < npc; j++)
                {
                    result += function(coefficients[i, 0], coefficients[j, 0]) * coefficients[i, 1] * coefficients[j, 1];
                }
            }
            return result ;
        }
    }
}
