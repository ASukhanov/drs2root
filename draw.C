{
TCanvas cc;
cc.Divide(2,2);
cc.cd(1);hch0->Draw("colz");
cc.cd(3);d4r_profile_0->Draw();
cc.cd(2);hch1->Draw("colz");
cc.cd(4);d4r_profile_1->Draw();
//cc.cd(2);tree->Draw("peak[1]");
//cc.cd(4);tree->Draw("peakpos[1]");
}
