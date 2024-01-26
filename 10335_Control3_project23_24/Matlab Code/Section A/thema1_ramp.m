  
    % Ramp Function
    opts = odeset('Refine',5);
    %Simulation 1
    
    [t,x] = ode45(@ramp_slope1,[0 5],[2; 1],opts);
    tc1 = tiledlayout(3,1);
    r = t;
    nexttile;
    plot(t,x(:,1),'-bo',t,x(:,2),'-r*')
    legend('x_1','x_2')
    ylabel('x_1,x_2');
    xlabel('t');
    title('Διάγραμμα Εξισώσεων κατάστασης')
    nexttile;
    x1 = x(:, 1);
    x2 = x(:, 2);
    plot(x1, x2);
    xlabel('x_1');
    ylabel('x_2');
    title('Phase Plane Plot');
    grid on;
    nexttile;
    y = r - x1;
    plot(t,r,'-bo',t,y,'-r*');
    xlabel('Time');
    legend('r','y')
    title('Απόκριση Εξόδου και Εισόδου Αναφοράς')
    title(tc1,'Συνάρτηση Ράμπας για το πρώτο σύστημα με αρχικές συνθήκες y(0)=-2,ydot(0) = 0')
    %============================================================================
    figure
    %Simulation 2
    tc1 = tiledlayout(3,1);
    [t,x] = ode45(@ramp_slope1,[0 5],[-1; 1],opts);
    r = t;
    nexttile;
    plot(t,x(:,1),'-bo',t,x(:,2),'-r*')
    legend('x_1','x_2')
    ylabel('x_1,x_2');
    xlabel('t');
    title('Διάγραμμα Εξισώσεων κατάστασης')
    nexttile;
    x1 = x(:, 1);
    x2 = x(:, 2);
    plot(x1, x2);
    xlabel('x_1');
    ylabel('x_2');
    title('Phase Plane Plot');
    grid on;
    nexttile;
    y = r - x1;
    plot(t,r,'-bo',t,y,'-r*');
    xlabel('Time');
    legend('r','y')
    title('Απόκριση Εξόδου και Εισόδου Αναφοράς')
    title(tc1,'Συνάρτηση Ράμπας για το πρώτο σύστημα με αρχικές συνθήκες y(0)=1,ydot(0) = 0')
    %============================================================================
    figure
    %Simulation 3
    tc1 = tiledlayout(3,1);
    [t,x] = ode45(@ramp_slope1,[0 5],[0; 0.5],opts);
    r = t;
    nexttile;
    plot(t,x(:,1),'-bo',t,x(:,2),'-r*')
    legend('x_1','x_2')
    ylabel('x_1,x_2');
    xlabel('t');
    title('Διάγραμμα Εξισώσεων κατάστασης')
    nexttile;
    x1 = x(:, 1);
    x2 = x(:, 2);
    plot(x1, x2);
    xlabel('x_1');
    ylabel('x_2');
    title('Phase Plane Plot');
    grid on;
    nexttile;
    y = r - x1;
    plot(t,r,'-bo',t,y,'-r*');
    xlabel('Time');
    legend('r','y')
    title('Απόκριση Εξόδου και Εισόδου Αναφοράς')
    title(tc1,'Συνάρτηση Ράμπας για το πρώτο σύστημα με αρχικές συνθήκες y(0)=0,ydot(0) = 0.5')
    %============================================================================
    figure
    %Simulation 4
    tc1 = tiledlayout(3,1);
    [t,x] = ode45(@ramp_slope1,[0 5],[-2; -1],opts);
    r = t;
    nexttile;
    plot(t,x(:,1),'-bo',t,x(:,2),'-r*')
    legend('x_1','x_2')
    ylabel('x_1,x_2');
    xlabel('t');
    title('Διάγραμμα Εξισώσεων κατάστασης')
    nexttile;
    x1 = x(:, 1);
    x2 = x(:, 2);
    plot(x1, x2);
    xlabel('x_1');
    ylabel('x_2');
    title('Phase Plane Plot');
    grid on;
    nexttile;
    y = r - x1;
    plot(t,r,'-bo',t,y,'-r*');
    xlabel('Time');
    legend('r','y')
    title('Απόκριση Εξόδου και Εισόδου Αναφοράς')
    title(tc1,'Συνάρτηση Ράμπας για το πρώτο σύστημα με αρχικές συνθήκες y(0)=2,ydot(0) = 2')
    %============================================================================
    figure
    %Simulation 5
    tc1 = tiledlayout(3,1);
    [t,x] = ode45(@ramp_slope1,[0 5],[-2.5; 2],opts);
    r = t;
    nexttile;
    plot(t,x(:,1),'-bo',t,x(:,2),'-r*')
    legend('x_1','x_2')
    ylabel('x_1,x_2');
    xlabel('t');
    title('Διάγραμμα Εξισώσεων κατάστασης')
    nexttile;
    x1 = x(:, 1);
    x2 = x(:, 2);
    plot(x1, x2);
    xlabel('x_1');
    ylabel('x_2');
    title('Phase Plane Plot');
    grid on;
    nexttile;
    y = r - x1;
    plot(t,r,'-bo',t,y,'-r*');
    xlabel('Time');
    legend('r','y')
    title('Απόκριση Εξόδου και Εισόδου Αναφοράς')
    title(tc1,'Συνάρτηση Ράμπας για το πρώτο σύστημα με αρχικές συνθήκες y(0)=2.5,ydot(0) = -1')
    %============================================================================
    figure
    %Simulation 6
    tc1 = tiledlayout(3,1);
    [t,x] = ode45(@ramp_slope1,[0 5],[-1.1; -1],opts);
    r = t;
    nexttile;
    plot(t,x(:,1),'-bo',t,x(:,2),'-r*')
    legend('x_1','x_2')
    ylabel('x_1,x_2');
    xlabel('t');
    title('Διάγραμμα Εξισώσεων κατάστασης')
    nexttile;
    x1 = x(:, 1);
    x2 = x(:, 2);
    plot(x1, x2);
    xlabel('x_1');
    ylabel('x_2');
    title('Phase Plane Plot');
    grid on;
    nexttile;
    y = r - x1;
    plot(t,r,'-bo',t,y,'-r*');
    xlabel('Time');
    legend('r','y')
    title('Απόκριση Εξόδου και Εισόδου Αναφοράς')
    title(tc1,'Συνάρτηση Ράμπας για το πρώτο σύστημα με αρχικές συνθήκες y(0)=1.1,ydot(0) = 2')
    %============================================================================

