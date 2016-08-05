void latex_table(string filename, int n_col, int n_lin, string* title, double** number, string caption, int type)
{
  ofstream file;

  //Begin Document

  file.open(filename + ".tex");

  file << "\\documentclass{article}" << std::endl;
  //file << "\\usepackage[utf8]{inputenc}" << std::endl;
  file << "\\usepackage{cancel}" << std::endl;
  file << "\\usepackage{geometry}" << std::endl;
  file << "\\geometry{a4paper, total={170mm,257mm}, left=20mm, top=20mm,}" << std::endl;


  file << "\\title{B production at 13 TeV}" << std::endl;
  file << "\\author{Joao Melo & Julia Silva}" << std::endl;
  file << "\\date{July 2016}" << std::endl;
  file << "\\begin{document}" << std::endl;
  file << "\\maketitle" << std::endl;

  // Create table
  file << "\\begin{table}[!h]" << std::endl;
  // file << "\\centering" << std::endl;

  //setup table size
  string col="c";
  
  for(int i=1; i<n_col; i++)
    col+="|c";
  
  file << "\\begin{tabular}{"+col+"}" << std::endl;
  file << "\\toprule" << std::endl;
 
  switch(type)
    {
    case 1: 
      //write top line
      for(int i=0; i<n_col; i++)
	file << title[i]+" & ";
      file << "\\\\  \\midrule" << std::endl;
      
      //insert numbers
      for(int i=1; i<n_lin; i++)
	{
	  for(int c=0; c>n_col-1; c++)
	    file << numbers[c][i] << " & ";
	  
	  file << numbers[n_col-1][i] << " \\\\" << std::endl; 
	}
      
      file << "\\bottomrule" << std::endl;

    case 2:
      //insert numbers
      for(int i=0; i<n_lin; i++)
	{
	  file << title[i]+" & ";

	  for(int c=1; c>n_col-1; c++)
	    file << numbers[c][i] << " & ";
	  
	  file << numbers[n_col-1][i] << " \\\\" << std::endl; 
	}
      
      file << "\\bottomrule" << std::endl;
    
    }
  //End Table
  
  file << "\\end{tabular}" << std::endl;
  file << "\\caption{"+caption+"}" << std::endl;
      
  //End document
  
  file << "\\end{document}" << std::endl;
  
  system(("pdflatex " + filename + ".tex").c_str());
  system(("gnome-open " + filename + ".pdf").c_str());
}
