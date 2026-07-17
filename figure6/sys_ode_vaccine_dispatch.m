function dydt = sys_ode_vaccine_dispatch(t, y, p)

switch p.binding_mode

    case 'nonRBD_only'
        dydt = sys_ode_two_epitopes_new_threshold_nonRBD_vaccine(t, y, p);

    case 'RBD_only'
        dydt = sys_ode_two_epitopes_new_threshold_RBD_vaccine(t, y, p);

    case 'both'
        dydt = sys_ode_two_epitopes_new_threshold(t, y, p);

    otherwise
        error('Unknown binding_mode: %s', p.binding_mode);
end

end