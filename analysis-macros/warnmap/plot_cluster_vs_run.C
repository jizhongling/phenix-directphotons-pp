plot_cluster_vs_run()
{

  TTree *tin = new TTree();
  tin->ReadFile("summaries/allruns_Run13pp510MinBias_cluster.txt","run/I:cluster/F");

  TCanvas *c1 = new TCanvas();
  tin->Draw("cluster");
  c1->Print("plots/cluster_per_event_hist.eps");
  c1->Print("plots/cluster_per_event_hist.png");

  TCanvas *c2 = new TCanvas();
  tin->Draw("cluster:Entry$");
  c2->Print("plots/cluster_per_event_graph.eps");
  c2->Print("plots/cluster_per_event_graph.png");

}
