function[numofmoves, success] = runtest(mapfile, robotstart, targetstart, payloads, fuels)

envmap = load(mapfile);

close all;

%draw the environment
figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(envmap); axis square; colorbar; colormap jet; hold on;

%current positions of the target and robot
robotpos = robotstart;
targetpos = targetstart;

% Constants
C = 10;
D = 10;
P = 20;
F = 50;
Fcur = 1000000;
Fmax = 1000000;

%now comes the main loop
hr = -1;
ht = -1;
hf = zeros(size(fuels,1));
hp = -1;
numofmoves = 0;
success = 0;

for i = 1:20000

    disp(size(fuels))
    disp(size(payloads))
    disp(size(fuels,1))

    %draw the positions
    if (hr ~= -1)
        delete(hr);
        delete(ht);
    end
    hr = text(robotpos(1), robotpos(2), 'R', 'Color', 'g', 'FontWeight', 'bold');
    ht = text(targetpos(1), targetpos(2), 'T', 'Color', 'm', 'FontWeight', 'bold');
    hr = scatter(robotpos(1), robotpos(2), 10, 'g', 'filled');
    ht = scatter(targetpos(1), targetpos(2), 10, 'm', 'filled');
    hp = text(payloads(1), payloads(2), 'P', 'Color', 'b', 'FontWeight', 'bold');
    for f = 1:size(fuels,1)
        hf(f) = text(fuels(f,1), fuels(f,2), 'F', 'Color', 'y', 'FontWeight', 'bold');
    end

    pause(0.1);
    %pause();

    %call robot planner to find what they want to do
    tStart = tic;
    action = robotplanner(envmap, robotpos, targetpos, payloads, fuels, C, D, P, F, Fcur, Fmax);

    %check that the new commanded position is valid
    if (newrobotpos(1) < 1 || newrobotpos(1) > size(envmap, 1) || ...
            newrobotpos(2) < 1 || newrobotpos(2) > size(envmap, 2))
        fprintf(1, 'ERROR: out-of-map robot position commanded\n');
        return;
    elseif (envmap(newrobotpos(1), newrobotpos(2)) ~= 0)
        fprintf(1, 'ERROR: invalid robot position commanded\n');
        return;
    end

    %make the moves

    actionD = cast(action,'like',robotpos);
    robotpos = robotpos + actionD;
    numofmoves = numofmoves + 1;
    fuelDelta = -C;
    if (actionD(0) ~= 0 || actionD(1) ~= 0)
        fuelDelta = fuelDelta - D;
        if(payloads(1) == robotpos(1) && payloads(2) == robotpos(2))
            fuelDelta = fuelDelta - P;
            payloads = payloads + actionD
        end
    else
        for f = 1:size(fuels,1)
            if(fuels(f,1) == robotpos(1) && fuels(f,2) == robotpos(2))
                fuelDelta = F;
            end
        end
    end

    % update fuel amount, check if fuel valid
    Fcur = Fcur + fuelDelta;
    if (Fcur > Fmax)
        Fcur = Fmax;
    end

    if (Fcur < 0)
        fprintf(1, 'ERROR: Robot ran out of fuel\n');
        return;
    end

    success = 1;
    %check if payloads are at target
    if (payloads(1) ~= targetpos(1) || payloads(2) ~= targetpos(2))
        success = 0;
    end

    if(success == 1)
        break;
    end
end

fprintf(1, 'success=%d number of moves made=%d\n', caught, numofmoves);
