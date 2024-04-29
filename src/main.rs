use cpal::{
    traits::{DeviceTrait, HostTrait, StreamTrait},
    FromSample, Sample, SizedSample,
};

fn main() {
    println!("Hello, world! {}...{}", u16::MIN, u16::MAX);

    let host = cpal::default_host();
    let device = host.default_output_device();
    if let Some(device) = device {    
        println!("Default host: {}", device.name().unwrap_or("null".to_string()));

        if let Ok(config) = device.default_output_config() {
            println!("SampleFormat: {}", config.sample_format());
            println!("SampleRate: {}", config.sample_rate().0);
            println!("Channels: {}", config.channels());
        }
    }
}
