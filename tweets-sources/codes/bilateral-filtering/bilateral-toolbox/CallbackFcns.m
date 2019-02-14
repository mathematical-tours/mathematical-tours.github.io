function CallbackFcns (action)

switch (action)
    case 'sigmas_slider'
        sigmas = get(gcbo, 'Value');
        sigmas = round(sigmas) + 1;
        updatesigmas(sigmas);
    case 'sigmar_slider'
        sigmar = get(gcbo, 'Value');
        sigmar = (round(sigmar*10))/10 + 5;
        updatesigmar(sigmar);
    case 'sigmas_edit'
        sigmas = eval(get(gcbo, 'String'));
        updatesigmas(sigmas);
    case 'sigmar_edit'
        sigmar = eval(get(gcbo, 'String'));
        updatesigmar(sigmar);
    case 'eps_slider'
        eps = get(gcbo, 'Value');
        eps = (round(eps) + 1)*0.001;
        update_eps(eps);
    case 'eps_edit'
        eps = eval(get(gcbo, 'String'));
        update_eps(eps);
    case 'browse'
        [filename, user_canceled] = imgetfile;
        if (~user_canceled)
            imagepath_Handle = findobj(gcbf,'Tag','imagepath');
            set(imagepath_Handle,'String',filename);
        end
    case 'load'
        imagepath_Handle = findobj(gcbf,'Tag','imagepath');
        mread = imread(get(imagepath_Handle,'String'));
        axes(findobj(gcbf,'Tag','input'));
        cla; hold on;
        imshow(mread); axis('image', 'off');
        minput = double(mread);
        set(gcbf,'UserData',minput);
    case 'filter'
        minput = get(gcbf,'UserData');
        editS_Handle = findobj(gcbf,'Tag','editS');
        sigmas = eval(get(editS_Handle,'String'));
        editR_Handle = findobj(gcbf,'Tag','editR');
        sigmar = eval(get(editR_Handle,'String'));
        editEps_Handle = findobj(gcbf,'Tag','editEps');
        eps = eval(get(editEps_Handle,'String'));
%         [moutput,params] = shiftableBF(minput,sigmas,sigmar);
        [moutput,N] = GPA(minput, sigmar, sigmas, eps, 'Gauss');
        axes(findobj(gcbf,'Tag','output'));
        cla; hold on;
        imshow(uint8(moutput)); axis('image', 'off');
        
        rkernel_Handle = findobj(gcbf,'Tag','rkernelplot');
        axes(rkernel_Handle);
        cla;
        set(rkernel_Handle,'Color','White');
        set(rkernel_Handle,'AmbientLightColor','White');
        set(rkernel_Handle,'XColor','Black');
        set(rkernel_Handle,'YColor','Black');
        box on;
        sigmar2 = sigmar^2;
        L = -127;
        U =128;
        t  = L : 0.01 : U;
        tauList = [-50,0,50];   %center
        for i=1:length(tauList)
            tau=tauList(i);
            g  = exp(-0.5*(t-tau).^2/sigmar2);
            % Gauss-polynomial
            nu = 0.5*(1/sigmar2);
            gapprox = zeros(size(t));
            for k = 0 : N
                gapprox = gapprox + (1/factorial(k)) * (1/sigmar2)^k ...
                    *(tau*t).^k;
            end
            gapprox = gapprox.*exp(-nu*t.^2).* exp(-nu*tau^2);
            %g3 = (1 - nu*((t-tau).^2/N)).^N;
            % display
            hold on, plot(t,g, 'r','LineWidth',3);
            hold on, plot(t,gapprox,'k','LineWidth',2);
            axis('tight'), grid('on'),
            hleg=legend('Target Gaussian', 'Gaussian-Polynomial');
            set(hleg,'fontsize',8);
        end
        
    case 'save'
        imsave(findobj(gcbf,'Tag','output'));
end

function updatesigmas (sigmas)
sliderS_Handle = findobj(gcbf,'Tag','sliderS');
set(sliderS_Handle,'Value',sigmas-1);
editS_Handle = findobj(gcbf,'Tag','editS');
set(editS_Handle,'String',sigmas);

function updatesigmar (sigmar)
sliderR_Handle = findobj(gcbf,'Tag','sliderR');
set(sliderR_Handle,'Value',sigmar-5);
editR_Handle = findobj(gcbf,'Tag','editR');
set(editR_Handle,'String',sigmar);

function update_eps (eps)
sliderEps_Handle = findobj(gcbf,'Tag','sliderEps');
set(sliderEps_Handle,'Value',eps*1000-1);
editEps_Handle = findobj(gcbf,'Tag','editEps');
set(editEps_Handle,'String',eps);