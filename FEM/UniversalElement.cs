namespace FEM
{
   public class UniversalElement
    {
       public double[,] integrationPoints = new double[4, 2];
      public double[,] N_Ksi = new double[4, 4];
      public double[,] N_Eta = new double[4, 4];
     public double[,] N = new double[4, 4];
       public double[,] IntegrationPoints2;

        public void PrepereMatrix()
        {
            integrationPoints[0, 0] = -0.577;
            integrationPoints[0, 1] = -0.577;
            integrationPoints[1, 0] = 0.577;
            integrationPoints[1, 1] = -0.577;
            integrationPoints[2, 0] = 0.577;
            integrationPoints[2, 1] = 0.577;
            integrationPoints[3, 0] = -0.577;
            integrationPoints[3, 1] = 0.577;

            IntegrationPoints2 = new double[8,2];
           
            //Dodawanie dodatkowych danych bo warunki brzegowe maja inne 1 i -1
            IntegrationPoints2[0,0] = -0.577;
            IntegrationPoints2[0,1] = -1;
            IntegrationPoints2[1,0] = 0.577;
            IntegrationPoints2[1,1] = -1;
            IntegrationPoints2[2,0] = 1;
            IntegrationPoints2[2,1] = -0.577;
            IntegrationPoints2[3,0] = 1;
            IntegrationPoints2[3,1] = 0.577;
            IntegrationPoints2[4,0] = 0.577;
            IntegrationPoints2[4,1] = 1;
            IntegrationPoints2[5,0] = -0.577;
            IntegrationPoints2[5,1] = 1;
            IntegrationPoints2[6,0] = -1;
            IntegrationPoints2[6,1] = 0.577;
            IntegrationPoints2[7,0] = -1;
            IntegrationPoints2[7,1] = -0.577;

            for (int i=0;i<4;i++)
            {
                int j = 0;
                N_Ksi[i, j] = -0.25 * (1 - integrationPoints[i, 1]);
                N_Ksi[i, j+1] = 0.25 * (1 - integrationPoints[i, 1]);
                N_Ksi[i, j+2] = 0.25 * (1 + integrationPoints[i, 1]);
                N_Ksi[i, j + 3] = -0.25 * (1 + integrationPoints[i, 1]);

                N_Eta[i,j]= -0.25 * (1 - integrationPoints[i, 0]);
                N_Eta[i, j+1] = -0.25 * (1 + integrationPoints[i, 0]);
                N_Eta[i, j+2]= 0.25 * (1 + integrationPoints[i, 0]);
                N_Eta[i, j+3]= 0.25 * (1 - integrationPoints[i, 0]);

                N[i, j] = 0.25 * (1 - integrationPoints[i, 0]) * (1 - integrationPoints[i, 1]);
                N[i, j+1] = 0.25 * (1 + integrationPoints[i, 0]) * (1 - integrationPoints[i, 1]);
                N[i, j+2] = 0.25 * (1 + integrationPoints[i, 0]) * (1 + integrationPoints[i, 1]);
                N[i, j+3] = 0.25 * (1 - integrationPoints[i, 0]) * (1 + integrationPoints[i, 1]);
            }
        }

    

        public double[,] CalculateJackobian(double[,] pointTab, int pointNr)
        {
            this.PrepereMatrix();
            double[,] resultTab = new double[2, 2];

            for(int i =0;i<4;i++)
            {
                resultTab[0, 0] += N_Ksi[pointNr, i] * pointTab[i, 0];
                resultTab[0, 1] += N_Ksi[pointNr, i] * pointTab[i, 1];
                resultTab[1, 0] += N_Eta[pointNr, i] * pointTab[i, 0];
                resultTab[1, 1] += N_Eta[pointNr, i] * pointTab[i, 1];
            }

            return resultTab;
        }
    }
}
