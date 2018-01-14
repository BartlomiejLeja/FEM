using System;
using System.Collections.Generic;

namespace FEM.FemGrid
{
    class Grid
    {
        public Grid()
        {
            GlobalData.GetGlobalData(@"E:\Studia\Informatyka_Stosowana\MES\cw\FEM\FEM\data.txt");
        }
     
        public List<Node> Nd = new List<Node>();
        public List<Element>  El = new  List<Element>();
        public List<double[,]> JakobianList = new List<double[,]>();
        const double eps = 1e-12; // stała przybliżenia zera
        public double[,] pointTab = new double[4, 2];
        double xPosition = 0.0;
        double yPosition = 0.0;
        int elementNumber=1;
        UniversalElement universalElement;

        public void FillFemGrid()
        {
            for(int i  = 0; i<GlobalData.ne;i++)
            {
                El.Add(new Element());
            }
            //Adding nodes
             universalElement = new UniversalElement();
            for (int i = 0; i < GlobalData.nB; i++)
            {
                for(int j = GlobalData.nH; j>0;j-- )
                {
                    if (i == 0 || i == GlobalData.nB  || j == 1 || j == GlobalData.nH )
                    {
                        
                        Nd.Add(new Node(xPosition, yPosition, true,100));//100 temperatura poczatkowa w kazdym wezle
                        yPosition += (GlobalData.H/(GlobalData.nH-1));
                    }
                    else
                    {
                        Nd.Add(new Node(xPosition, yPosition, false,100));
                        yPosition += (GlobalData.H / (GlobalData.nH-1 ));
                    }
                  
                }
                yPosition = 0.0;
                xPosition += (GlobalData.B/(GlobalData.nB-1));
            }
            //adding elements
            //for (int i = 0; i< GlobalData.nB -1; i++)
            //{
            //    for (int j = GlobalData.nH-1; j > 0; j--)
            //    {
            //        if (j == 1)
            //        {
            //            El.Add(new Element(elementNumber, elementNumber + GlobalData.nH, elementNumber + GlobalData.nH + 1, elementNumber + 1));
            //            elementNumber += 2;
            //        }
            //        else
            //        {
            //            El.Add(new Element(elementNumber, elementNumber + GlobalData.nH , elementNumber + GlobalData.nH + 1, elementNumber + 1));
            //            elementNumber += 1;
            //        }
            //    }
            //}

            int m = 1;
            for (int i = 1, j = 1; i <= GlobalData.ne; i++, j++)
            {
               
                El[i - 1].IdTab[0] = j;
                El[i - 1].IdTab[3] = j + 1;
                if (j == m * GlobalData.nH - 1)
                {
                    j = m * GlobalData.nH;
                    m++;
                }
            }
            m = 2;
            for (int i = 1, j = GlobalData.nH + 1; i <= GlobalData.ne; i++, j++)
            {
                El[i - 1].IdTab[1] = j;
                El[i - 1].IdTab[2] = j + 1;
                if (j == m * GlobalData.nH - 1)
                {
                    j = m * GlobalData.nH;
                    m++;
                }
            }

            for (int i = 0; i < GlobalData.ne; i++) //nn zamiast ne?
            {
                Nd[El[i].IdTab[0] - 1].Temp = universalElement.N[0,0] * Nd[El[i].IdTab[0] - 1].Temp + universalElement.N[0,1] * Nd[El[i].IdTab[1] - 1].Temp + universalElement.N[0,2] * Nd[El[i].IdTab[2] - 1].Temp + universalElement.N[0,3] * Nd[El[i].IdTab[3] - 1].Temp;
                Nd[El[i].IdTab[1] - 1].Temp = universalElement.N[1,0] * Nd[El[i].IdTab[0] - 1].Temp + universalElement.N[1,1] * Nd[El[i].IdTab[1] - 1].Temp + universalElement.N[1,2] * Nd[El[i].IdTab[2] - 1].Temp + universalElement.N[1,3] * Nd[El[i].IdTab[3] - 1].Temp;
                Nd[El[i].IdTab[2] - 1].Temp = universalElement.N[2,0] * Nd[El[i].IdTab[0] - 1].Temp + universalElement.N[2,1] * Nd[El[i].IdTab[1] - 1].Temp + universalElement.N[2,2] * Nd[El[i].IdTab[2] - 1].Temp + universalElement.N[2,3] * Nd[El[i].IdTab[3] - 1].Temp;
                Nd[El[i].IdTab[3] - 1].Temp = universalElement.N[3,0] * Nd[El[i].IdTab[0] - 1].Temp + universalElement.N[3,1] * Nd[El[i].IdTab[1] - 1].Temp + universalElement.N[3,2] * Nd[El[i].IdTab[2] - 1].Temp + universalElement.N[3,3] * Nd[El[i].IdTab[3] - 1].Temp;
            }

            for (int j = 0; j < El.Count; j++)
            {
                for (int i = 0; i < 4; i++) // dla każdego elemntu 4 jakobiany dla 4 punktów całkowania
                {
                    pointTab[i, 0] = Nd[El[j].IdTab[i]-1].X;
                    pointTab[i, 1] = Nd[El[j].IdTab[i]-1].Y;
                }
                for(int k=0;k<4;k++) // 36 jakobianów bo 9 razy 4
                {
                    //4 jakobiany bo 4 punkty calkowania tablica [2,2] bo jabobian 2d
                    JakobianList.Add( universalElement.CalculateJackobian(pointTab, k));
                }
            }
            
        }

        double [,] calculate_H1_matrix(double[,] jacobiMatrix, double detJ, UniversalElement Element, int pointNr) //obliczam pierwszą część macierzy H
        {
            double[] dNdx = new double[4]; //strona 2
            double[] dNdy = new double[4];
            double x = 0, y = 0;

            double[,] temp = new double[4,4];

            for (int i = 0; i < 4; i++)
            {
                dNdx[i] = (1 / detJ) * (jacobiMatrix[0,0] * Element.N_Ksi[pointNr,i] + jacobiMatrix[0,1] * Element.N_Eta[pointNr,i]); //strona pierwsza 
                dNdy[i] = (1 / detJ) * (jacobiMatrix[1,0] * Element.N_Ksi[pointNr,i] + jacobiMatrix[1,1] * Element.N_Eta[pointNr,i]);
            }

            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    x = dNdx[i] * dNdx[j];
                    y = dNdy[i] * dNdy[j];
                    temp[i,j] = (x + y) * GlobalData.k * detJ;//mnożymy po detJ bo po objętości przekstalcenie z lokalnego do globalnego
                }
            }
            return temp;
        }

        double[,] calculate_H2_matrix(double p1, double p2, double L) //obliczam druga część macierzy H
        {
            double[] N = new double[4];

            double[,] H2 = new double[4,4];
           
            N[0] = 0.25 * (1 - p1) * (1 - p2);
            N[1] = 0.25 * (1 + p1) * (1 - p2);
            N[2] = 0.25 * (1 + p1) * (1 + p2);
            N[3] = 0.25 * (1 - p1) * (1 + p2);

            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {

                    H2[i,j] = (N[i] * N[j]) * GlobalData.alfa * (L / 2);
                }
            }
            return H2;
        }

        double[] calculate_P_vector(double p1, double p2, double L) //obliczam wektor P
        {
            double[] N = new double[4];
            double[] P = new double[4];

            N[0] = 0.25 * (1 - p1) * (1 - p2);
            N[1] = 0.25 * (1 + p1) * (1 - p2);
            N[2] = 0.25 * (1 + p1) * (1 + p2);
            N[3] = 0.25 * (1 - p1) * (1 + p2);

            for (int i = 0; i < 4; i++)
            {
                P[i] = -GlobalData.alfa * -GlobalData.temp * N[i] * (L / 2);
            }

            return P;
        }

        double[,] calculate_C_matrix(double[,] jacobiMatrix, double detJ, UniversalElement Element, int pointNr) //obliczam macierz C po 4 punktach calkowania
        { 
            double[,] C = new double[4,4];
            
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    C[i,j] = (Element.N[pointNr,i] * Element.N[pointNr,j]) * GlobalData.c * GlobalData.ro * detJ;
                }
            }
            return C;
        }

        bool gauss(int n, double[,] AB, double[] X)//bo czasmi sie moze nie dac jak sa na przekatnych zera X wektor rozwiazan
        {
            int i, j, k;
            double m, s;

            for (i = 0; i < n - 1; i++)
            {
                for (j = i + 1; j < n; j++)
                {
                    if (Math.Abs(AB[i,i]) < eps) return false;
                    m = -AB[j,i] / AB[i,i];
                    for (k = i + 1; k <= n; k++)
                        AB[j,k] += m * AB[i,k];
                }
            }

            for (i = n - 1; i >= 0; i--)
            {
                s = AB[i,n];
                for (j = n - 1; j >= i + 1; j--)
                    s -= AB[i,j] * X[j];
                if (Math.Abs(AB[i,i]) < eps) return false;
                X[i] = s / AB[i,i];
            }
            return true;
        }

        double[,] matrixPlusMatrix(double[,] A, double[,] B, int size)
        {
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                {
                    A[i,j] += B[i,j];
                }
            return A;
        }

        double[] vectorPlusVector(double[] A, double[] B, int size)
        {
            for (int i = 0; i < size; i++)
                A[i] += B[i];

            return A;
        }

        public void commuteEveryThing()
        {
            double[,] p = new double[4,2];
           
            double[] detJ = new double[GlobalData.ne];

            double[,] tempH1 = new double[4,4];
           
            double[,] tempH2 = new double[4,4];
            
            double[] tempP1 = new double[4];

            double[] tempP2 = new double[4];
           
            double[] PG = new double[GlobalData.nH * GlobalData.nB];
           
            double[,] HG = new double[GlobalData.nH * GlobalData.nB, GlobalData.nB * GlobalData.nH];

            double[,] CG = new double[GlobalData.nH * GlobalData.nB, GlobalData.nB * GlobalData.nH];//test

            double[,] CL = new double[4, 4];//test

            double[,] GaussMatrix = new double[GlobalData.nH * GlobalData.nB, GlobalData.nB * GlobalData.nH + 1];
           
            double[,] HL = new double[4,4]; //HL to zssumowane wszystkie h z jednego punktu calkowania
           
            int deltaT = 50; //takie dane

            double[,] H2 = new double[4,4];
           
            double[] P = new double[4];

            double[] PL = new double[4];

            double[] Ct = new double[4];

            double[,] result = new double[2,2];
            
            double[,] H1 = new double[4,4];
           
            double[,] C = new double[4,4];
            
            double[] t0 = new double[4];
           
            double[] t1 = new double[16];


           
            //for (int i = 0; i <= 500; i+=deltaT)  //petla po czasie 
            //{
            for (int j = 0; j < GlobalData.ne; j++) //petla po elementach
            {
                for (int i=0;i<4;i++)
                {
                    for(int s=0;s<4;s++)
                    {
                        H1[i, s] = 0;
                        H2[i, s] = 0;
                        HL[i, s] = 0;
                        CL[i, s] = 0;
                        C[i, s] = 0;
                    }
                    PL[i] = 0;
                    Ct[i] = 0;
                  //  t0[i] = 0;
                }
                
                p[0,0] = Nd[El[j].IdTab[0] - 1].X;
                p[0,1] = Nd[El[j].IdTab[0] - 1].Y;
                p[1,0] = Nd[El[j].IdTab[1] - 1].X;
                p[1,1] = Nd[El[j].IdTab[1] - 1].Y;
                p[2,0] = Nd[El[j].IdTab[2] - 1].X;
                p[2,1] = Nd[El[j].IdTab[2] - 1].Y;
                p[3,0] = Nd[El[j].IdTab[3] - 1].X;
                p[3,1] = Nd[El[j].IdTab[3] - 1].Y;

                for (int l = 0; l < 4; l++)
                {
                    t0[l] = Nd[El[j].IdTab[l] - 1].Temp;
                }

                for (int k = 0; k < 4; k++) //pętla po 4 punktach całkowania
                {
                    for (int i = 0; i < 4; i++)
                    {
                        Ct[i] = 0;
                    }
                    result = JakobianList[j*4+k];
                    detJ[j] = result[0,0] * result[1,1] - result[0,1] * result[1,0];
                    H1 = calculate_H1_matrix(result, detJ[j], universalElement, k);
                    C = calculate_C_matrix(result, detJ[j], universalElement, k);
                    for (int i = 0; i < 4; i++)
                    {
                        for (int t = 0; t < 4; t++)
                        {
                            C[i, t] /= deltaT;
                        }
                    }

                    for (int i = 0; i < 4; i++)
                    {
                        for (int z = 0; z < 4; z++)
                        {
                            Ct[i] += C[i,z] * t0[z];
                        }
                    }

                    PL = vectorPlusVector(PL, Ct, 4);// P z daszkiem

                    HL = matrixPlusMatrix(HL, H1, 4);
                   // CL = matrixPlusMatrix(CL, C, 4);
                   HL = matrixPlusMatrix(HL, C, 4);//H z daszkiem

                }

                //sprawdzanie 4 ifami czy jest warunek brzegowy moze byc X=0 
                //petla po 2 punktach calkowania warunki brzegowe
                if (Nd[El[j].IdTab[0] - 1].Y == 0 || Nd[El[j].IdTab[1] - 1].Y == 0)
                {
                    //liczenie calki po powieszchni czyli 1D wyznacznik macierzy dx/2 2 bo w ukladzie lokalnym dlugosc jest od -1 do 1
                    for (int i = 0; i < 4; i++)
                    {
                        for (int l = 0; l < 4; l++)
                        {
                            tempH1[i, l] = 0;
                            tempH2[i, l] = 0;
                            H2[i, l] = 0;
                        }
                        P[i] = 0;
                        tempP1[i] = 0;
                        tempP2[i] = 0;
                    }


                    tempH1 = calculate_H2_matrix(universalElement.IntegrationPoints2[0,0], universalElement.IntegrationPoints2[0,1], GlobalData.L1);
                    H2 = matrixPlusMatrix(H2, tempH1, 4);
                    tempH2 = calculate_H2_matrix(universalElement.IntegrationPoints2[1,0], universalElement.IntegrationPoints2[1,1], GlobalData.L1);
                    H2 = matrixPlusMatrix(H2, tempH1, 4);
                    HL = matrixPlusMatrix(HL, H2, 4);
                    tempP1 = calculate_P_vector(universalElement.IntegrationPoints2[0,0], universalElement.IntegrationPoints2[0,1], GlobalData.L1);
                    P = vectorPlusVector(P, tempP1, 4);
                    tempP2 = calculate_P_vector(universalElement.IntegrationPoints2[1,0], universalElement.IntegrationPoints2[1,1], GlobalData.L1);
                    P = vectorPlusVector(P, tempP2, 4);
                    PL = vectorPlusVector(PL,P, 4);
                }
                if (Nd[El[j].IdTab[1] - 1].X == GlobalData.B || Nd[El[j].IdTab[2] - 1].X == GlobalData.B)
                {
                    for (int i = 0; i < 4; i++)
                    {
                        for (int l = 0; l < 4; l++)
                        {
                            tempH1[i, l] = 0;
                            tempH2[i, l] = 0;
                            H2[i, l] = 0;
                        }
                        P[i] = 0;
                        tempP1[i] = 0;
                        tempP2[i] = 0;
                    }

                    tempH1 = calculate_H2_matrix(universalElement.IntegrationPoints2[2,0], universalElement.IntegrationPoints2[2,1], GlobalData.L2);
                    H2 = matrixPlusMatrix(H2, tempH1, 4);
                    tempH2 = calculate_H2_matrix(universalElement.IntegrationPoints2[3,0], universalElement.IntegrationPoints2[3,1], GlobalData.L2);
                    H2 = matrixPlusMatrix(H2, tempH1, 4);
                    HL = matrixPlusMatrix(HL, H2, 4);
                    tempP1 = calculate_P_vector(universalElement.IntegrationPoints2[2,0], universalElement.IntegrationPoints2[2,1], GlobalData.L2);
                    P = vectorPlusVector(P, tempP1, 4);
                    tempP2 = calculate_P_vector(universalElement.IntegrationPoints2[3,0], universalElement.IntegrationPoints2[3,1], GlobalData.L2);
                    P = vectorPlusVector(P, tempP2, 4);
                    PL = vectorPlusVector(PL, P, 4);
                }

                if (Nd[El[j].IdTab[2] - 1].Y == GlobalData.H || Nd[El[j].IdTab[3] - 1].Y == GlobalData.H)
                {
                    for (int i = 0; i < 4; i++)
                    {
                        for (int l = 0; l < 4; l++)
                        {
                            tempH1[i, l] = 0;
                            tempH2[i, l] = 0;
                            H2[i, l] = 0;
                        }
                        P[i] = 0;
                        tempP1[i] = 0;
                        tempP2[i] = 0;
                    }

                    tempH1 = calculate_H2_matrix(universalElement.IntegrationPoints2[4,0], universalElement.IntegrationPoints2[4,1], GlobalData.L1);
                    H2 = matrixPlusMatrix(H2, tempH1, 4);
                    tempH2 = calculate_H2_matrix(universalElement.IntegrationPoints2[5,0], universalElement.IntegrationPoints2[5,1], GlobalData.L1);
                    H2 = matrixPlusMatrix(H2, tempH1, 4);
                    HL = matrixPlusMatrix(HL, H2, 4);
                    tempP1 = calculate_P_vector(universalElement.IntegrationPoints2[4,0], universalElement.IntegrationPoints2[4,1], GlobalData.L1);
                    P = vectorPlusVector(P, tempP1, 4);
                    tempP2 = calculate_P_vector(universalElement.IntegrationPoints2[5,0], universalElement.IntegrationPoints2[5,1], GlobalData.L1);
                    P = vectorPlusVector(P, tempP2, 4);
                    PL = vectorPlusVector(PL, P, 4);
                }
                if (Nd[El[j].IdTab[3] - 1].X == 0 || Nd[El[j].IdTab[0] - 1].X == 0)
                {
                    for (int i = 0; i < 4; i++)
                    {
                        for (int l = 0; l < 4; l++)
                        {
                            tempH1[i, l] = 0;
                            tempH2[i, l] = 0;
                            H2[i, l] = 0;
                        }
                        P[i] = 0;
                        tempP1[i] = 0;
                        tempP2[i] = 0;
                    }

                    tempH1 = calculate_H2_matrix(universalElement.IntegrationPoints2[6,0], universalElement.IntegrationPoints2[6,1], GlobalData.L2);
                    H2 = matrixPlusMatrix(H2, tempH1, 4);
                    tempH2 = calculate_H2_matrix(universalElement.IntegrationPoints2[7,0], universalElement.IntegrationPoints2[7,1], GlobalData.L2);
                    H2 = matrixPlusMatrix(H2, tempH1, 4);
                    HL = matrixPlusMatrix(HL, H2, 4);
                    tempP1 = calculate_P_vector(universalElement.IntegrationPoints2[6,0], universalElement.IntegrationPoints2[6,1], GlobalData.L2);
                    P = vectorPlusVector(P, tempP1, 4);
                    tempP2 = calculate_P_vector(universalElement.IntegrationPoints2[7,0], universalElement.IntegrationPoints2[7,1], GlobalData.L2);
                    P = vectorPlusVector(P, tempP2, 4);
                    PL = vectorPlusVector(PL, P, 4);
                }

                for (int i = 0; i < 4; i++)
                {
                    for (int k = 0; k < 4; k++)
                    {
                        int l = El[j].IdTab[i] - 1;
                        int m = El[j].IdTab[k] - 1;
                        HG[l,m] += HL[i,k];// agregacja
                        CG[El[j].IdTab[i] - 1, El[j].IdTab[k] - 1] +=CL[i, k];
                    }
                }

                for (int i = 0; i < 4; i++)
                {
                    PG[El[j].IdTab[i] - 1] += PL[i];
                }
            } //koniec pętli po elementach

            for (int i = 0; i < 16; i++)
            {
                for (int j = 0; j < 16; j++)
                {
                    GaussMatrix[i,j] = HG[i,j];
                }
            }
            for (int k = 0; k < 16; k++)
            {
                GaussMatrix[k,16] = PG[k];
            }

            gauss(16, GaussMatrix, t1);
        }
    }
}
