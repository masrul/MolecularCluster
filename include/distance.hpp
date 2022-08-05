/**
 * Author      : Masrul Huda (mail2masrul@gmail.com)
 * Host        : iMacPro@Swalm
 * Created     : Friday Aug 05, 2022 16:01:52 CDT
 */

void min_img_dr(float* dr,matrix& box){
    for (int k=0;k<3;++k){
        while (dr[k]>0.5*box[k][k] || dr[k]<-0.5*box[k][k]){
            if (dr[k] > 0.5*box[k][k])
                dr[k] -=box[k][k];
            else if (dr[k] < -0.5*box[k][k])
                dr[k] +=box[k][k];
        }
    }
}
