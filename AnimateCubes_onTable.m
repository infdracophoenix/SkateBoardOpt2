function [] =AnimateCubes_onTable(XList,XcList,OList,OcList,F_Mat,Fc_Mat,Ftable_Mat,P_base,Pc_base,Ptable_base,S,frameRate,fileName,rotate,axlims)

NumFrames=size(XList,4);
NumCubes=size(XList,3);

figHandle=figure();
    for j=1:NumCubes
        
        Verts=(S*0.5*OList(:,:,j,1)*P_base+XList(:,:,j,1))';
        
       H(j)= patch('Vertices',Verts,'Faces',F_Mat,'FaceVertexCData',hsv(6),'FaceColor','flat');
        axis equal
        axis([-.3 .3 -.3 .3 0 .1])
        view([-37.5 30])
    end
        Verts=(OcList(:,:,1)*Pc_base+XcList(:,:,1))';
        H(NumCubes+1)= patch('Vertices',Verts,'Faces',Fc_Mat,'FaceVertexCData',hsv(2),'FaceColor','flat','FaceVertexAlphaData',[0;1],'FaceAlpha','flat');
        
        Verts=Ptable_base;
        H(NumCubes+2)= patch('Vertices',Verts,'Faces',Ftable_Mat,'FaceVertexCData',[1 1 1],'FaceColor','flat');
        
    viewAzim=linspace(-37.5,-37.5-180,NumFrames);
for k=1:NumFrames
    Verts=(OcList(:,:,k)*Pc_base+XcList(:,:,k))';
    H(NumCubes+1).Vertices=Verts;
    
    for j=1:NumCubes
        
        Verts=(S*0.5*OList(:,:,j,k)*P_base+XList(:,:,j,k))';
        H(j).Vertices=Verts;
        axis equal
        axis(axlims)
        if rotate
        view([viewAzim(k) 30])
        end
        grid on
        
    end

    pause(1/frameRate);
    drawnow
    FlipBook(k)=getframe(figHandle,[10,10,520,400]);
    
end

Writer=VideoWriter(fileName,'MPEG-4');

Writer.FrameRate = frameRate;

open(Writer);
writeVideo(Writer,FlipBook);
close(Writer);





end