namespace FEM.FemGrid
{
    class Element
    {
        public int[] IdTab = new int[4];
        public double[] H;

        public Element()
        {
          
        }
        public Element(int Id1,int Id2, int Id3, int Id4)
        {
            IdTab[0] = Id1;
            IdTab[1] = Id2;
            IdTab[2] = Id3;
            IdTab[3] = Id4;
        }

        public void GenerateHArray()
        {

        }
    }
}
