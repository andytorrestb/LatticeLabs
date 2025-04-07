clc;
clear all;

nframes = 50;
L = 2 * pi; % Define domain length
[x, y] = meshgrid(linspace(0, L, 20), linspace(0, L, 20)); % Grid for quiver

v = VideoWriter('quiver_animation.avi');
open(v)

figure
ax = gca();
axis([0 L 0 L])
axis equal

for t = 1:nframes
    % Define vector field (example: rotating flow)
    U = cos(x) .* sin(y + 0.1 * t);
    V = -sin(x + 0.1 * t) .* cos(y);
    
    if t == 1
        h = quiver(ax, x, y, U, V, 'LineWidth', 2);
        set(ax, 'XLimMode', 'manual', 'YLimMode', 'manual');
    else
        set(h, 'UData', U, 'VData', V);
    end
    
    drawnow();
    % writeVideo(v, getframe(ax));
end

close(v);

disp(size(x));
disp(size(U));