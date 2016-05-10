{
TCanvas cc;
cc.Divide(2,2);
cc.cd(1);hch1->Draw("colz");
cc.cd(3);d4r_profile_1->Draw();
cc.cd(2);tree->Draw("peak[1]");
cc.cd(4);tree->Draw("peakpos[1]");
}
